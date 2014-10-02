//------------------------------------------------------------------------------
// Includes
//------------------------------------------------------------------------------
#include "MantidAlgorithms/MonteCarloAbsorption.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/WorkspaceValidators.h"
#include "MantidAPI/SampleEnvironment.h"
#include "MantidGeometry/IDetector.h"
#include "MantidGeometry/Objects/Track.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/NeutronAtom.h"
#include "MantidKernel/VectorHelper.h"

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

/// @cond
namespace
{
  /// Number of attempts to choose a random point within the object before it gives up
  const int MaxRandPointAttempts(500);

  /// Element size in mm
  const double ELEMENT_SIZE = 1.0;
}
/// @endcond

namespace Mantid
{
  namespace Algorithms
  {

    DECLARE_ALGORITHM(MonteCarloAbsorption)

    using namespace API;
    using namespace Geometry;
    using namespace Kernel;

    //------------------------------------------------------------------------------
    // Public methods
    //------------------------------------------------------------------------------

    /**
     * Constructor
     */
    MonteCarloAbsorption::MonteCarloAbsorption() :
      m_samplePos(), m_sourcePos(), m_blocks(),
      m_blkHalfX(0.0), m_blkHalfY(0.0), m_blkHalfZ(0.0),
      m_rng(NULL), m_inputWS(), m_sampleShape(NULL), m_container(NULL), m_numberOfPoints(0),
      m_xStepSize(0), m_numberOfEvents(300)
    {
    }

    /**
     * Destructor
     */
    MonteCarloAbsorption::~MonteCarloAbsorption()
    {
      delete m_rng;
    }

    //------------------------------------------------------------------------------
    // Private methods
    //------------------------------------------------------------------------------

    /**
     * Initialize the algorithm
     */
    void MonteCarloAbsorption::init()
    {
      // The input workspace must have an instrument and units of wavelength
      auto wsValidator = boost::make_shared<CompositeValidator>();
      wsValidator->add<WorkspaceUnitValidator>("Wavelength");
      wsValidator->add<InstrumentValidator>();

      declareProperty(new WorkspaceProperty<>("InputWorkspace", "", Direction::Input,
                                              wsValidator),
                      "The name of the input workspace.  The input workspace must have X units of wavelength.");
      declareProperty(new WorkspaceProperty<> ("OutputWorkspace", "", Direction::Output),
                      "The name to use for the output workspace.");
      auto positiveInt = boost::make_shared<Kernel::BoundedValidator<int> >();
      positiveInt->setLower(1);
      declareProperty("NumberOfWavelengthPoints", EMPTY_INT(), positiveInt,
                      "The number of wavelength points for which a simulation is atttempted (default: all points)");
      declareProperty("EventsPerPoint", m_numberOfEvents, positiveInt,
                      "The number of \"neutron\" events to generate per simulated point");
      declareProperty("SeedValue", 123456789, positiveInt,
                      "Seed the random number generator with this value");

    }

    /**
     * Execution code
     */
    void MonteCarloAbsorption::exec()
    {
      retrieveInput();
      initCaches();

      g_log.debug() << "Creating output workspace\n";
      MatrixWorkspace_sptr correctionFactors = WorkspaceFactory::Instance().create(m_inputWS);
      correctionFactors->isDistribution(true); // The output of this is a distribution
      correctionFactors->setYUnit(""); // Need to explicitly set YUnit to nothing
      correctionFactors->setYUnitLabel("Attenuation factor");

      const bool isHistogram = m_inputWS->isHistogramData();
      const int numHists = static_cast<int>(m_inputWS->getNumberHistograms());
      const int numBins = static_cast<int>(m_inputWS->blocksize());

      // Compute the step size
      m_xStepSize = numBins/m_numberOfPoints;

      g_log.information() << "Simulation performed every " << m_xStepSize
                          << " wavelength points" << std::endl;

      Progress prog(this,0.0,1.0, numHists*numBins/m_xStepSize);
      for( int i = 0; i < numHists; ++i )
      {
        // Copy over the X-values
        const MantidVec & xValues = m_inputWS->readX(i);
        correctionFactors->dataX(i) = xValues;
        // Final detector position
        IDetector_const_sptr detector;
        try
        {
          detector = m_inputWS->getDetector(i);
        }
        catch(Kernel::Exception::NotFoundError&)
        {
          continue;
        }

        MantidVec & yValues = correctionFactors->dataY(i);
        MantidVec & eValues = correctionFactors->dataE(i);
        // Simulation for each requested wavelength point
        for( int bin = 0; bin < numBins; bin += m_xStepSize )
        {
          prog.report("Computing corrections for bin " + boost::lexical_cast<std::string>(bin));
          const double lambda = isHistogram ?
                (0.5 * (xValues[bin] + xValues[bin + 1]) ) : xValues[bin];
          doSimulation(detector.get(), lambda, yValues[bin], eValues[bin]);
          // Ensure we have the last point for the interpolation
          if ( m_xStepSize > 1 && bin + m_xStepSize >= numBins && bin+1 != numBins)
          {
            bin = numBins - m_xStepSize - 1;
          }
        }

        // Interpolate through points not simulated
        if( m_xStepSize > 1 )
        {
          prog.report("Interpolating unsimulated points");
          Kernel::VectorHelper::linearlyInterpolateY(xValues, yValues, m_xStepSize);
        }

      }

      // Save the results
      setProperty("OutputWorkspace", correctionFactors);
    }


    /**
     * Perform the simulation
     * @param detector :: A pointer to the current detector
     * @param lambda :: The chosen wavelength
     * @param attenFactor :: [Output] The calculated attenuation factor for this wavelength
     * @param error :: [Output] The value of the error on the factor
     */
    void MonteCarloAbsorption::doSimulation(const IDetector *const detector, const double lambda,
                                            double & attenFactor, double & error)
    {
      /**
       Currently, assuming square beam profile to pick start position then randomly selecting
       a point within the sample using it's bounding box.
       This point defines the single scattering point and hence the attenuation path lengths and final
       directional vector to the detector
       */
      // Absolute detector position
      const V3D detectorPos(detector->getPos());

      int numDetected(0);
      attenFactor = 0.0;
      error = 0.0;
      while( numDetected < m_numberOfEvents )
      {
        V3D startPos = sampleBeamProfile();
        V3D scatterPoint = selectScatterPoint();
        try
        {
          attenFactor += attenuationFactor(startPos, scatterPoint, detectorPos, lambda);
        }
        catch(std::logic_error &)
        {
          continue;
        }
        ++numDetected;
      }

      // Attenuation factor is simply the average value
      attenFactor /= numDetected;
      // Error is 1/sqrt(nevents)
      error = 1./sqrt((double)numDetected);
    }

    /**
     * Sample the beam profile for a random start location, assigning a weight
     * to the given point
     * @returns The starting location
     */
    V3D MonteCarloAbsorption::sampleBeamProfile() const
    {
      return m_sourcePos;
    }

    /**
     * Selects a random location within the sample + container environment. The bounding box is
     * used as an approximation to generate a point and this is then tested for its validity within
     * the shape.
     * @returns Selected position as V3D object
     */
    V3D MonteCarloAbsorption::selectScatterPoint() const
    {
      // Randomly select a block from the subdivided set and then randomly select a point
      // within that block and test if it inside the sample/container. If yes then accept, else
      // keep trying.
      boost::uniform_int<size_t> uniIntDist(0, m_blocks.size() - 1);
      boost::variate_generator<boost::mt19937&, boost::uniform_int<size_t>> uniInt(*m_rng, uniIntDist);
      boost::uniform_real<> uniRealDist(0, 1.0);
      boost::variate_generator<boost::mt19937&, boost::uniform_real<>> uniReal(*m_rng, uniRealDist);

      V3D scatterPoint;
      int nattempts(0);
      while( nattempts < MaxRandPointAttempts )
      {
        const auto & block = m_blocks[uniInt()];
        const double x = m_blkHalfX*(2.0*uniReal() - 1.0) + block.xMin();
        const double y = m_blkHalfY*(2.0*uniReal() - 1.0) + block.yMin();
        const double z = m_blkHalfZ*(2.0*uniReal() - 1.0) + block.zMin();
        scatterPoint(x,y,z);
        ++nattempts;
        if( ptIntersectsSample(scatterPoint) )
        {
          scatterPoint += m_samplePos;
          return scatterPoint;
        }
      }
      // If we got here then the shape is too strange for the bounding box to be of any use.
      g_log.error() << "Attempts to generate a random point with the sample/can "
                    << "have exceeded the allowed number of tries.\n";
      throw std::runtime_error("Attempts to produce random scatter point failed. Check sample shape.");
    }

    /**
     * Return the attenuation factor for the given track
     * @param startPos :: The origin of the track
     * @param scatterPoint :: The point of scatter
     * @param finalPos :: The end point of the track
     * @param lambda :: The wavelength of the neutron
     * @returns The attenuation factor for this neutron's track
     */
    double
    MonteCarloAbsorption::attenuationFactor(const V3D & startPos, const V3D & scatterPoint,
                                            const V3D & finalPos, const double lambda)
    {
      double factor(1.0);

      // Define two tracks, before and after scatter, and trace check their
      // intersections with the the environment and sample
      Track beforeScatter(scatterPoint, (startPos - scatterPoint));
      Track afterScatter(scatterPoint, (finalPos - scatterPoint));
      // Theoretically this should never happen as there should always be an intersection
      // but do to precision limitations points very close to the surface give
      // zero intersection, so just reject
      if( m_sampleShape->interceptSurface(beforeScatter) == 0 ||
          m_sampleShape->interceptSurface(afterScatter) == 0 )
      {
        throw std::logic_error("Track has no surfaces intersections.");
      }

      double length = beforeScatter.begin()->distInsideObject;
      factor *= attenuation(length, *m_sampleMaterial, lambda);

      beforeScatter.clearIntersectionResults();
      if( m_container )
      {
        m_container->interceptSurfaces(beforeScatter);
      }
      // Attenuation factor is product of factor for each material
      Track::LType::const_iterator cend = beforeScatter.end();
      for(Track::LType::const_iterator citr = beforeScatter.begin();
          citr != cend; ++citr)
      {
        length = citr->distInsideObject;
        IObjComponent *objComp = dynamic_cast<IObjComponent*>(citr->componentID);
        Material_const_sptr mat = objComp->material();
        factor *= attenuation(length, *mat, lambda);
      }

      length = afterScatter.begin()->distInsideObject;
      factor *= attenuation(length, *m_sampleMaterial, lambda);

      afterScatter.clearIntersectionResults();
      if( m_container )
      {
        m_container->interceptSurfaces(afterScatter);
      }
      // Attenuation factor is product of factor for each material
      cend = afterScatter.end();
      for(Track::LType::const_iterator citr = afterScatter.begin();
          citr != cend; ++citr)
      {
        length = citr->distInsideObject;
        IObjComponent *objComp = dynamic_cast<IObjComponent*>(citr->componentID);
        Material_const_sptr mat = objComp->material();
        factor *= attenuation(length, *mat, lambda);
      }

      return factor;
    }

    /**
     * Calculate the attenuation for a given length, material and wavelength
     * @param length :: Distance through the material
     * @param material :: A reference to the Material
     * @param lambda :: The wavelength
     * @returns The attenuation factor
     */
    double
    MonteCarloAbsorption::attenuation(const double length, const Kernel::Material& material,
                                      const double lambda) const
    {
      const double rho = material.numberDensity() * 100.0;
      const double sigma_s = material.totalScatterXSection(lambda);
      const double sigma_t = sigma_s + material.absorbXSection(lambda);
      return exp(-rho*sigma_t*length);
    }

    /**
     * Gather the input values and check validity
     * @throw std::invalid_argument If the input is invalid. Currently if there is
     * no defined sample shape
     */
    void MonteCarloAbsorption::retrieveInput()
    {
      m_inputWS = getProperty("InputWorkspace");

      m_sampleShape = &(m_inputWS->sample().getShape());
      m_sampleMaterial = &(m_inputWS->sample().getMaterial());
      if( !m_sampleShape->hasValidShape() )
      {
        g_log.debug() << "Invalid shape defined on workspace. TopRule = " << m_sampleShape->topRule()
                      << ", No. of surfaces: " << m_sampleShape->getSurfacePtr().size() << "\n";
        throw std::invalid_argument("Input workspace has an invalid sample shape.");
      }

      if( m_sampleMaterial->totalScatterXSection(1.0) == 0.0 )
      {
        g_log.warning() << "The sample material appears to have zero scattering cross section.\n"
                        << "Result will most likely be nonsensical.\n";
      }

      try
      {
        m_container = &(m_inputWS->sample().getEnvironment());
      }
      catch(std::runtime_error&)
      {
        m_container = NULL;
        g_log.information() << "No environment has been defined, continuing with only sample.\n";
      }

      m_numberOfPoints = getProperty("NumberOfWavelengthPoints");
      if( isEmpty(m_numberOfPoints) ||  m_numberOfPoints > static_cast<int>(m_inputWS->blocksize()) )
      {
        m_numberOfPoints = static_cast<int>(m_inputWS->blocksize());
        if( !isEmpty(m_numberOfPoints) )
        {
          g_log.warning() << "The requested number of wavelength points is larger than the spectra size. "
                          << "Defaulting to spectra size.\n";
        }
      }

      m_numberOfEvents = getProperty("EventsPerPoint");
    }

    /**
     * Initialise the caches used here including setting up the random
     * number generator
     */
    void MonteCarloAbsorption::initCaches()
    {
      g_log.debug() << "Caching input\n";
      if( m_rng ) delete m_rng;
      const int seedValue = getProperty("SeedValue");
      m_rng = new boost::mt19937(seedValue);

      m_samplePos = m_inputWS->getInstrument()->getSample()->getPos();
      m_sourcePos = m_inputWS->getInstrument()->getSource()->getPos();
      BoundingBox box(m_sampleShape->getBoundingBox());
      if( m_container )
      {
        BoundingBox envBox;
        m_container->getBoundingBox(envBox);
        box.grow(envBox);
      }

      // Chop the bounding box up into a set of small boxes. This will be used
      // as a first guess for generating a random scatter point
      const double cubeSide = ELEMENT_SIZE*1e-3;
      const double xLength = box.width().X();
      const double yLength = box.width().Y();
      const double zLength = box.width().Z();
      const int numXSlices = static_cast<int>(xLength/cubeSide);
      const int numYSlices = static_cast<int>(yLength/cubeSide);
      const int numZSlices = static_cast<int>(zLength/cubeSide);
      const double xThick = xLength/numXSlices;
      const double yThick = yLength/numYSlices;
      const double zThick = zLength/numZSlices;
      m_blkHalfX = 0.5*xThick;
      m_blkHalfY = 0.5*yThick;
      m_blkHalfZ = 0.5*zThick;

      const size_t numVolumeElements = numXSlices*numYSlices*numZSlices;
      g_log.debug() << "Attempting to divide sample + container into " << numVolumeElements << " blocks.\n";

      try
      {
        m_blocks.clear();
        m_blocks.reserve(numVolumeElements);
      }
      catch (std::exception&)
      {
        // Typically get here if the number of volume elements is too large
        // Provide a bit more information
        g_log.error("Too many volume elements requested - try increasing the value of the ElementSize property.");
        throw;
      }
      const auto boxCentre = box.centrePoint();
      const double x0 = boxCentre.X() - 0.5*xLength;
      const double y0 = boxCentre.Y() - 0.5*yLength;
      const double z0 = boxCentre.Z() - 0.5*zLength;
      // Store a chunk as a BoundingBox object.
      // Only cache blocks that have some intersection with the
      // sample or container.

      for (int i = 0; i < numZSlices; ++i)
      {
        const double zmin = z0 + i*zThick;
        const double zmax = zmin + zThick;
        for (int j = 0; j < numYSlices; ++j)
        {
          const double ymin = y0 + j*yThick;
          const double ymax = ymin + yThick;
          for (int k = 0; k < numXSlices; ++k)
          {
            const double xmin = x0 + k*xThick;
            const double xmax = xmin + xThick;
            BoundingBox block(xmax, ymax, zmax, xmin, ymin, zmin);
            if(boxIntersectsSample(block)) m_blocks.push_back(block);
          }
        }
      }

      g_log.debug() << "Sample + container divided into " << m_blocks.size() << " blocks.";
      if(m_blocks.size() == numVolumeElements) g_log.debug("\n");
      else g_log.debug() << " Skipped " << (numVolumeElements-m_blocks.size())
                         << " blocks that do not intersect with the sample + container\n";
    }

    /**
     * @param block A BoundingBox object defining a block
     * @return True if any of the vertices intersect the sample or container
     */
    bool MonteCarloAbsorption::boxIntersectsSample(const BoundingBox &block) const
    {
      // Check all 8 corners for intersection
      V3D vertices[] = {
        V3D(block.xMax(), block.yMin(), block.zMin()), // left-front-bottom
        V3D(block.xMax(), block.yMax(), block.zMin()), // left-front-top
        V3D(block.xMin(), block.yMax(), block.zMin()), // right-front-top
        V3D(block.xMin(), block.yMin(), block.zMin()), // right-front-bottom
        V3D(block.xMax(), block.yMin(), block.zMax()), // left-back-bottom
        V3D(block.xMax(), block.yMax(), block.zMax()), // left-back-top
        V3D(block.xMin(), block.yMax(), block.zMax()), // right-back-top
        V3D(block.xMin(), block.yMin(), block.zMax()) // right-back-bottom
      };
      bool intersects(true);
      for(size_t i = 0; i < 8; ++i)
      {
        intersects &= ptIntersectsSample(vertices[i]);
      }
      return intersects;
    }

    /**
     *
     * @param pt A V3D giving a point to test
     * @return True if point is inside sample or container, false otherwise
     */
    bool MonteCarloAbsorption::ptIntersectsSample(const V3D &pt) const
    {
      return m_sampleShape->isValid(pt) ||
          (m_container && m_container->isValid(pt));
    }

  }
}
