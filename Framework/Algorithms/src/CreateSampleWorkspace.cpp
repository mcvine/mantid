#include "MantidAlgorithms/CreateSampleWorkspace.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/FunctionDomain1D.h"
#include "MantidAPI/FunctionProperty.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidGeometry/Objects/ShapeFactory.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidGeometry/Instrument/RectangularDetector.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidKernel/MersenneTwister.h"

#include <cmath>
#include <ctime>
#include <numeric>
#include <stdexcept>

namespace Mantid {
namespace Algorithms {
using namespace Kernel;
using namespace API;
using namespace Geometry;
using namespace DataObjects;
using Mantid::MantidVec;
using Mantid::MantidVecPtr;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(CreateSampleWorkspace)

//----------------------------------------------------------------------------------------------
/** Constructor
 */
CreateSampleWorkspace::CreateSampleWorkspace() : m_randGen(nullptr) {}

//----------------------------------------------------------------------------------------------
/** Destructor
 */
CreateSampleWorkspace::~CreateSampleWorkspace() { delete m_randGen; }

//----------------------------------------------------------------------------------------------
/// Algorithm's name for identification. @see Algorithm::name
const std::string CreateSampleWorkspace::name() const {
  return "CreateSampleWorkspace";
}

/// Algorithm's version for identification. @see Algorithm::version
int CreateSampleWorkspace::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string CreateSampleWorkspace::category() const {
  return "Utility\\Workspaces";
}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void CreateSampleWorkspace::init() {
  declareProperty(make_unique<WorkspaceProperty<>>("OutputWorkspace", "",
                                                   Direction::Output),
                  "An output workspace.");
  std::vector<std::string> typeOptions{"Histogram", "Event"};
  declareProperty("WorkspaceType", "Histogram",
                  boost::make_shared<StringListValidator>(typeOptions),
                  "The type of workspace to create (default: Histogram)");

  // pre-defined function strings these use $PCx$ to define peak centres values
  // that will be replaced before use
  //$PC0$ is the far left of the data, and $PC10$ is the far right, and
  // therefore will often not be used
  //$PC5$ is the centre of the data
  m_preDefinedFunctionmap.emplace(
      "One Peak", "name=LinearBackground, A0=0.3; name=Gaussian, "
                  "PeakCentre=$PC5$, Height=10, Sigma=0.7;");
  m_preDefinedFunctionmap.emplace(
      "Multiple Peaks", "name=LinearBackground, A0=0.3;name=Gaussian, "
                        "PeakCentre=$PC3$, Height=10, Sigma=0.7;name=Gaussian, "
                        "PeakCentre=$PC6$, Height=8, Sigma=0.5");
  m_preDefinedFunctionmap.emplace("Flat background",
                                  "name=LinearBackground, A0=1;");
  m_preDefinedFunctionmap.emplace("Exp Decay",
                                  "name=ExpDecay, Height=100, Lifetime=1000;");
  m_preDefinedFunctionmap.emplace(
      "Powder Diffraction",
      "name= LinearBackground,A0=0.0850208,A1=-4.89583e-06;"
      "name=Gaussian,Height=0.584528,PeakCentre=$PC1$,Sigma=14.3772;"
      "name=Gaussian,Height=1.33361,PeakCentre=$PC2$,Sigma=15.2516;"
      "name=Gaussian,Height=1.74691,PeakCentre=$PC3$,Sigma=15.8395;"
      "name=Gaussian,Height=0.950388,PeakCentre=$PC4$,Sigma=19.8408;"
      "name=Gaussian,Height=1.92185,PeakCentre=$PC5$,Sigma=18.0844;"
      "name=Gaussian,Height=3.64069,PeakCentre=$PC6$,Sigma=19.2404;"
      "name=Gaussian,Height=2.8998,PeakCentre=$PC7$,Sigma=21.1127;"
      "name=Gaussian,Height=2.05237,PeakCentre=$PC8$,Sigma=21.9932;"
      "name=Gaussian,Height=8.40976,PeakCentre=$PC9$,Sigma=25.2751;");
  m_preDefinedFunctionmap.emplace(
      "Quasielastic", "name=Lorentzian,FWHM=0.3,PeakCentre=$PC5$,Amplitude=0.8;"
                      "name=Lorentzian,FWHM=0.1,PeakCentre=$PC5$,Amplitude=1;"
                      "name=LinearBackground,A0=0.1");
  m_preDefinedFunctionmap.emplace(
      "Quasielastic Tunnelling",
      "name=LinearBackground,A0=0.1;"
      "name=Lorentzian,FWHM=0.1,PeakCentre=$PC5$,Amplitude=1;"
      "name=Lorentzian,FWHM=0.05,PeakCentre=$PC7$,Amplitude=0.04;"
      "name=Lorentzian,FWHM=0.05,PeakCentre=$PC3$,Amplitude=0.04;"
      "name=Lorentzian,FWHM=0.05,PeakCentre=$PC8$,Amplitude=0.02;"
      "name=Lorentzian,FWHM=0.05,PeakCentre=$PC2$,Amplitude=0.02");
  m_preDefinedFunctionmap.emplace("User Defined", "");
  std::vector<std::string> functionOptions;
  functionOptions.reserve(m_preDefinedFunctionmap.size());
  for (const auto &preDefinedFunction : m_preDefinedFunctionmap) {
    functionOptions.push_back(preDefinedFunction.first);
  }
  declareProperty("Function", "One Peak",
                  boost::make_shared<StringListValidator>(functionOptions),
                  "Preset options of the data to fill the workspace with");
  declareProperty(
      "UserDefinedFunction", "",
      "Parameters defining the fitting function and its initial values");

  declareProperty("NumBanks", 2,
                  boost::make_shared<BoundedValidator<int>>(0, 100),
                  "The Number of banks in the instrument (default:2)");
  declareProperty("BankPixelWidth", 10,
                  boost::make_shared<BoundedValidator<int>>(0, 10000),
                  "The number of pixels in horizontally and vertically in a "
                  "bank (default:10)");
  declareProperty("NumEvents", 1000,
                  boost::make_shared<BoundedValidator<int>>(0, 100000),
                  "The number of events per detector, this is only used for "
                  "EventWorkspaces (default:1000)");
  declareProperty(
      "Random", false,
      "Whether to randomise the placement of events and data (default:false)");

  declareProperty("XUnit", "TOF",
                  "The unit to assign to the XAxis (default:\"TOF\")");
  declareProperty("XMin", 0.0, "The minimum X axis value (default:0)");
  declareProperty("XMax", 20000.0, "The maximum X axis value (default:20000)");
  declareProperty("BinWidth", 200.0,
                  boost::make_shared<BoundedValidator<double>>(0, 100000, true),
                  "The bin width of the X axis (default:200)");
  declareProperty("PixelSpacing", 0.008,
                  boost::make_shared<BoundedValidator<double>>(0, 100000, true),
                  "The spacing between detector pixels in M (default:0.008)");
  declareProperty("BankDistanceFromSample", 5.0,
                  boost::make_shared<BoundedValidator<double>>(0, 1000, true),
                  "The distance along the beam direction from the sample to "
                  "bank in M (default:5.0)");
  declareProperty("SourceDistanceFromSample", 10.0,
                  boost::make_shared<BoundedValidator<double>>(0, 1000, true),
                  "The distance along the beam direction from the source to "
                  "the sample in M (default:10.0)");
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void CreateSampleWorkspace::exec() {
  const std::string wsType = getProperty("WorkspaceType");
  const std::string preDefinedFunction = getProperty("Function");
  const std::string userDefinedFunction = getProperty("UserDefinedFunction");
  const int numBanks = getProperty("NumBanks");
  const int bankPixelWidth = getProperty("BankPixelWidth");
  const int numEvents = getProperty("NumEvents");
  const bool isRandom = getProperty("Random");
  const std::string xUnit = getProperty("XUnit");
  const double xMin = getProperty("XMin");
  const double xMax = getProperty("XMax");
  double binWidth = getProperty("BinWidth");
  const double pixelSpacing = getProperty("PixelSpacing");
  const double bankDistanceFromSample = getProperty("BankDistanceFromSample");
  const double sourceSampleDistance = getProperty("SourceDistanceFromSample");

  if (xMax <= xMin) {
    throw std::invalid_argument("XMax must be larger than XMin");
  }

  if (binWidth > (xMax - xMin)) {
    // the bin width is so large that there is less than one bin - so adjust it
    // down
    binWidth = xMax - xMin;
    g_log.warning() << "The bin width is so large that there is less than one "
                       "bin - it has been changed to " << binWidth << std::endl;
  }

  std::string functionString = "";
  if (m_preDefinedFunctionmap.find(preDefinedFunction) !=
      m_preDefinedFunctionmap.end()) {
    // extract pre-defined string
    functionString = m_preDefinedFunctionmap[preDefinedFunction];
  }
  if (functionString.empty()) {
    functionString = userDefinedFunction;
  }

  if (!m_randGen) {
    int seedValue = 0;
    if (isRandom) {
      seedValue = static_cast<int>(std::time(nullptr));
    }
    m_randGen = new Kernel::MersenneTwister(seedValue);
  }

  int numPixels = numBanks * bankPixelWidth * bankPixelWidth;

  Progress progress(this, 0, 1, numBanks);

  // Create an instrument with one or more rectangular banks.
  Instrument_sptr inst = createTestInstrumentRectangular(
      progress, numBanks, bankPixelWidth, pixelSpacing, bankDistanceFromSample,
      sourceSampleDistance);

  int num_bins = static_cast<int>((xMax - xMin) / binWidth);

  MatrixWorkspace_sptr ws;
  if (wsType == "Event") {
    ws = createEventWorkspace(numPixels, num_bins, numEvents, xMin, binWidth,
                              bankPixelWidth * bankPixelWidth, inst,
                              functionString, isRandom);
  } else {
    ws = createHistogramWorkspace(numPixels, num_bins, xMin, binWidth,
                                  bankPixelWidth * bankPixelWidth, inst,
                                  functionString, isRandom);
  }
  // add chopper
  this->addChopperParameters(ws);

  // Set the Unit of the X Axis
  try {
    ws->getAxis(0)->unit() = UnitFactory::Instance().create(xUnit);
  } catch (Exception::NotFoundError &) {
    ws->getAxis(0)->unit() = UnitFactory::Instance().create("Label");
    Unit_sptr unit = ws->getAxis(0)->unit();
    boost::shared_ptr<Units::Label> label =
        boost::dynamic_pointer_cast<Units::Label>(unit);
    label->setLabel(xUnit, xUnit);
  }

  ws->setYUnit("Counts");
  ws->setTitle("Test Workspace");
  DateAndTime run_start("2010-01-01T00:00:00");
  DateAndTime run_end("2010-01-01T01:00:00");
  Run &theRun = ws->mutableRun();
  // belt and braces use both approaches for setting start and end times
  theRun.setStartAndEndTime(run_start, run_end);
  theRun.addLogData(new PropertyWithValue<std::string>(
      "run_start", run_start.toISO8601String()));
  theRun.addLogData(
      new PropertyWithValue<std::string>("run_end", run_end.toISO8601String()));

  // Assign it to the output workspace property
  setProperty("OutputWorkspace", ws);
  ;
}
/** Add chopper to the existing matrix workspace
@param ws  -- shared pointer to existing matrix workspace which has instrument
and chopper

@returns workspace modified to have Fermi chopper added to it.
*/
void CreateSampleWorkspace::addChopperParameters(
    API::MatrixWorkspace_sptr &ws) {

  auto testInst = ws->getInstrument();
  auto chopper = testInst->getComponentByName("chopper-position");

  // add chopper parameters
  auto &paramMap = ws->instrumentParameters();
  const std::string description(
      "The initial rotation phase of the disk used to calculate the time"
      " for neutrons arriving at the chopper according to the formula time = "
      "delay + initial_phase/Speed");
  paramMap.add<double>("double", chopper.get(), "initial_phase", -3000.,
                       &description);
  paramMap.add<std::string>("string", chopper.get(), "ChopperDelayLog",
                            "fermi_delay");
  paramMap.add<std::string>("string", chopper.get(), "ChopperSpeedLog",
                            "fermi_speed");
  paramMap.add<std::string>("string", chopper.get(), "FilterBaseLog",
                            "is_running");
  paramMap.add<bool>("bool", chopper.get(), "filter_with_derivative", false);
}

/** Create histogram workspace
 */
MatrixWorkspace_sptr CreateSampleWorkspace::createHistogramWorkspace(
    int numPixels, int numBins, double x0, double binDelta,
    int start_at_pixelID, Geometry::Instrument_sptr inst,
    const std::string &functionString, bool isRandom) {
  MantidVecPtr x, y, e;
  x.access().resize(numBins + 1);
  e.access().resize(numBins);
  for (int i = 0; i < numBins + 1; ++i) {
    x.access()[i] = x0 + i * binDelta;
  }

  std::vector<double> xValues(x.access().begin(), x.access().end() - 1);
  y.access() = evalFunction(functionString, xValues, isRandom ? 1 : 0);
  e.access().resize(numBins);

  // calculate e as sqrt(y)
  typedef double (*uf)(double);
  uf dblSqrt = std::sqrt;
  std::transform(y.access().begin(), y.access().end(), e.access().begin(),
                 dblSqrt);

  MatrixWorkspace_sptr retVal = boost::make_shared<DataObjects::Workspace2D>();
  retVal->initialize(numPixels, numBins + 1, numBins);
  retVal->setInstrument(inst);

  for (size_t wi = 0; wi < static_cast<size_t>(numPixels); wi++) {
    retVal->setX(wi, x);
    retVal->setData(wi, y, e);
    retVal->getSpectrum(wi)->setDetectorID(detid_t(start_at_pixelID + wi));
    retVal->getSpectrum(wi)->setSpectrumNo(specnum_t(wi + 1));
  }

  return retVal;
}

/** Create event workspace
 */
EventWorkspace_sptr CreateSampleWorkspace::createEventWorkspace(
    int numPixels, int numBins, int numEvents, double x0, double binDelta,
    int start_at_pixelID, Geometry::Instrument_sptr inst,
    const std::string &functionString, bool isRandom) {
  DateAndTime run_start("2010-01-01T00:00:00");

  // add one to the number of bins as this is histogram
  int numXBins = numBins + 1;

  auto retVal = boost::make_shared<EventWorkspace>();
  retVal->initialize(numPixels, 1, 1);

  retVal->setInstrument(inst);

  // Create the x-axis for histogramming.
  MantidVecPtr x1;
  MantidVec &xRef = x1.access();
  xRef.resize(numXBins);
  for (int i = 0; i < numXBins; ++i) {
    xRef[i] = x0 + i * binDelta;
  }

  // Set all the histograms at once.
  retVal->setAllX(x1);

  std::vector<double> xValues(xRef.begin(), xRef.end() - 1);
  std::vector<double> yValues =
      evalFunction(functionString, xValues, isRandom ? 1 : 0);

  // we need to normalise the results and then multiply by the number of events
  // to find the events per bin
  double sum_of_elems = std::accumulate(yValues.begin(), yValues.end(), 0.0);
  double event_distrib_factor = numEvents / sum_of_elems;
  std::transform(yValues.begin(), yValues.end(), yValues.begin(),
                 std::bind1st(std::multiplies<double>(), event_distrib_factor));
  // the array should now contain the number of events required per bin

  // Make fake events
  size_t workspaceIndex = 0;

  const double hourInSeconds = 60 * 60;
  for (int wi = 0; wi < numPixels; wi++) {
    EventList &el = retVal->getEventList(workspaceIndex);
    el.setSpectrumNo(wi + 1);
    el.setDetectorID(wi + start_at_pixelID);

    // for each bin

    for (int i = 0; i < numBins; ++i) {
      // create randomised events within the bin to match the number required -
      // calculated in yValues earlier
      int eventsInBin = static_cast<int>(yValues[i]);
      for (int q = 0; q < eventsInBin; q++) {
        DateAndTime pulseTime =
            run_start + (m_randGen->nextValue() * hourInSeconds);
        el += TofEvent((i + m_randGen->nextValue()) * binDelta, pulseTime);
      }
    }
    workspaceIndex++;
  }

  return retVal;
}
//----------------------------------------------------------------------------------------------
/**
 * Evaluates a function and returns the values as a vector
 *
 *
 * @param functionString :: the function string
 * @param xVal :: A vector of the x values
 * @param noiseScale :: A scaling factor for niose to be added to the data, 0=
 *no noise
 * @returns the calculated values
 */
std::vector<double>
CreateSampleWorkspace::evalFunction(const std::string &functionString,
                                    const std::vector<double> &xVal,
                                    double noiseScale = 0) {
  size_t xSize = xVal.size();
  // replace $PCx$ values
  std::string parsedFuncString = functionString;
  for (int x = 0; x <= 10; ++x) {
    // get the rough peak centre value
    int index = static_cast<int>((xSize / 10) * x);
    if ((x == 10) && (index > 0))
      --index;
    double replace_val = xVal[index];
    std::ostringstream tokenStream;
    tokenStream << "$PC" << x << "$";
    std::string token = tokenStream.str();
    std::string replaceStr = boost::lexical_cast<std::string>(replace_val);
    replaceAll(parsedFuncString, token, replaceStr);
  }

  IFunction_sptr func_sptr =
      FunctionFactory::Instance().createInitialized(parsedFuncString);
  FunctionDomain1DVector fd(xVal);
  FunctionValues fv(fd);
  func_sptr->function(fd, fv);

  std::vector<double> results;
  results.resize(xSize);
  for (size_t x = 0; x < xSize; ++x) {
    results[x] = fv.getCalculated(x);
    if (noiseScale != 0) {
      results[x] += ((m_randGen->nextValue() - 0.5) * noiseScale);
    }
    // no negative values please  - it messes up the error calculation
    results[x] = fabs(results[x]);
  }
  return results;
}

void CreateSampleWorkspace::replaceAll(std::string &str,
                                       const std::string &from,
                                       const std::string &to) {
  if (from.empty())
    return;
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing
                              // 'x' with 'yx'
  }
}

//----------------------------------------------------------------------------------------------
/**
 * Create an test instrument with n panels of rectangular detectors,
 *pixels*pixels in size,
 * a source and spherical sample shape.
 *
 * Banks' lower-left corner is at position (0,0,5*banknum) and they go up to
 *(pixels*0.008, pixels*0.008, Z)
 * Pixels are 4 mm wide.
 *
 * @param progress :: progress indicator
 * @param num_banks :: number of rectangular banks to create
 * @param pixels :: number of pixels in each direction.
 * @param pixelSpacing :: padding between pixels
 * @param bankDistanceFromSample :: Distance of first bank from sample (defaults
 *to 5.0m)
 * @param sourceSampleDistance :: The distance from the source to the sample
 * @returns A shared pointer to the generated instrument
 */
Instrument_sptr CreateSampleWorkspace::createTestInstrumentRectangular(
    API::Progress &progress, int num_banks, int pixels, double pixelSpacing,
    const double bankDistanceFromSample, const double sourceSampleDistance) {
  auto testInst = boost::make_shared<Instrument>("basic_rect");
  // The instrument is going to be set up with z as the beam axis and y as the
  // vertical axis.
  testInst->setReferenceFrame(
      boost::make_shared<ReferenceFrame>(Y, Z, Left, ""));

  const double cylRadius(pixelSpacing / 2);
  const double cylHeight(0.0002);
  // One object
  Object_sptr pixelShape = createCappedCylinder(
      cylRadius, cylHeight, V3D(0.0, -cylHeight / 2.0, 0.0), V3D(0., 1.0, 0.),
      "pixel-shape");

  for (int banknum = 1; banknum <= num_banks; banknum++) {
    // Make a new bank
    std::ostringstream bankname;
    bankname << "bank" << banknum;

    RectangularDetector *bank = new RectangularDetector(bankname.str());
    bank->initialize(pixelShape, pixels, 0.0, pixelSpacing, pixels, 0.0,
                     pixelSpacing, banknum * pixels * pixels, true, pixels);

    // Mark them all as detectors
    for (int x = 0; x < pixels; x++)
      for (int y = 0; y < pixels; y++) {
        boost::shared_ptr<Detector> detector = bank->getAtXY(x, y);
        if (detector)
          // Mark it as a detector (add to the instrument cache)
          testInst->markAsDetector(detector.get());
      }

    testInst->add(bank);
    // Set the bank along the z-axis of the instrument. (beam direction).
    bank->setPos(V3D(0.0, 0.0, bankDistanceFromSample * banknum));

    progress.report();
  }

  // Define a source component
  ObjComponent *source =
      new ObjComponent("moderator", Object_sptr(new Object), testInst.get());
  source->setPos(V3D(0.0, 0.0, -sourceSampleDistance));
  testInst->add(source);
  testInst->markAsSource(source);

  // Add chopper
  ObjComponent *chopper = new ObjComponent(
      "chopper-position", Object_sptr(new Object), testInst.get());
  chopper->setPos(V3D(0.0, 0.0, -0.25 * sourceSampleDistance));
  testInst->add(chopper);

  // Define a sample as a simple sphere
  Object_sptr sampleSphere =
      createSphere(0.001, V3D(0.0, 0.0, 0.0), "sample-shape");
  ObjComponent *sample =
      new ObjComponent("sample", sampleSphere, testInst.get());
  testInst->setPos(0.0, 0.0, 0.0);
  testInst->add(sample);
  testInst->markAsSamplePos(sample);

  return testInst;
}

//----------------------------------------------------------------------------------------------
/**
 * Create a capped cylinder object
 */
Object_sptr CreateSampleWorkspace::createCappedCylinder(double radius,
                                                        double height,
                                                        const V3D &baseCentre,
                                                        const V3D &axis,
                                                        const std::string &id) {
  std::ostringstream xml;
  xml << "<cylinder id=\"" << id << "\">"
      << "<centre-of-bottom-base x=\"" << baseCentre.X() << "\" y=\""
      << baseCentre.Y() << "\" z=\"" << baseCentre.Z() << "\"/>"
      << "<axis x=\"" << axis.X() << "\" y=\"" << axis.Y() << "\" z=\""
      << axis.Z() << "\"/>"
      << "<radius val=\"" << radius << "\" />"
      << "<height val=\"" << height << "\" />"
      << "</cylinder>";

  ShapeFactory shapeMaker;
  return shapeMaker.createShape(xml.str());
}

//----------------------------------------------------------------------------------------------

/**
 * Create a sphere object
 */
Object_sptr CreateSampleWorkspace::createSphere(double radius,
                                                const V3D &centre,
                                                const std::string &id) {
  ShapeFactory shapeMaker;
  std::ostringstream xml;
  xml << "<sphere id=\"" << id << "\">"
      << "<centre x=\"" << centre.X() << "\"  y=\"" << centre.Y() << "\" z=\""
      << centre.Z() << "\" />"
      << "<radius val=\"" << radius << "\" />"
      << "</sphere>";
  return shapeMaker.createShape(xml.str());
}

} // namespace Algorithms
} // namespace Mantid
