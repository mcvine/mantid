//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidCurveFitting/Algorithms/Fit.h"
#include "MantidCurveFitting/CostFunctions/CostFuncFitting.h"

#include "MantidAPI/CostFunctionFactory.h"
#include "MantidAPI/FuncMinimizerFactory.h"
#include "MantidAPI/FunctionValues.h"
#include "MantidAPI/IFuncMinimizer.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/TableRow.h"
#include "MantidAPI/WorkspaceFactory.h"

#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/StartsWithValidator.h"

//=====================================================
#include "MantidCurveFitting/GSLVector.h"
#include "MantidCurveFitting/FuncMinimizers/LocalSearchMinimizer.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidKernel/ArrayProperty.h"

namespace Mantid {
namespace CurveFitting {
namespace Algorithms {

// Register the class into the algorithm factory
DECLARE_ALGORITHM(Fit)

/** Initialisation method
*/
void Fit::initConcrete() {

  declareProperty("Ties", "", Kernel::Direction::Input);
  getPointerToProperty("Ties")
      ->setDocumentation("Math expressions defining ties between parameters of "
                         "the fitting function.");
  declareProperty("Constraints", "", Kernel::Direction::Input);
  getPointerToProperty("Constraints")->setDocumentation("List of constraints");
  auto mustBePositive = boost::shared_ptr<Kernel::BoundedValidator<int>>(
      new Kernel::BoundedValidator<int>());
  mustBePositive->setLower(0);
  declareProperty(
      "MaxIterations", 500, mustBePositive->clone(),
      "Stop after this number of iterations if a good fit is not found");
  declareProperty("OutputStatus", "", Kernel::Direction::Output);
  getPointerToProperty("OutputStatus")
      ->setDocumentation("Whether the fit was successful");
  declareProperty("OutputChi2overDoF", 0.0, "Returns the goodness of the fit",
                  Kernel::Direction::Output);

  // Disable default gsl error handler (which is to call abort!)
  gsl_set_error_handler_off();

  std::vector<std::string> minimizerOptions =
      API::FuncMinimizerFactory::Instance().getKeys();

  declareProperty("Minimizer", "Levenberg-Marquardt",
                  Kernel::IValidator_sptr(
                      new Kernel::StartsWithValidator(minimizerOptions)),
                  "Minimizer to use for fitting. Minimizers available are "
                  "\"Levenberg-Marquardt\", \"Simplex\", \"FABADA\", "
                  "\"Conjugate gradient (Fletcher-Reeves imp.)\", \"Conjugate "
                  "gradient (Polak-Ribiere imp.)\", \"BFGS\", and "
                  "\"Levenberg-MarquardtMD\"");

  std::vector<std::string> costFuncOptions =
      API::CostFunctionFactory::Instance().getKeys();
  // select only CostFuncFitting variety
  for (auto it = costFuncOptions.begin(); it != costFuncOptions.end(); ++it) {
    auto costFunc = boost::dynamic_pointer_cast<CostFunctions::CostFuncFitting>(
        API::CostFunctionFactory::Instance().create(*it));
    if (!costFunc) {
      *it = "";
    }
  }
  declareProperty(
      "CostFunction", "Least squares",
      Kernel::IValidator_sptr(
          new Kernel::ListValidator<std::string>(costFuncOptions)),
      "The cost function to be used for the fit, default is Least squares",
      Kernel::Direction::InOut);
  declareProperty(
      "CreateOutput", false,
      "Set to true to create output workspaces with the results of the fit"
      "(default is false).");
  declareProperty(
      "Output", "",
      "A base name for the output workspaces (if not "
      "given default names will be created). The "
      "default is to use the name of the original data workspace as prefix "
      "followed by suffixes _Workspace, _Parameters, etc.");
  declareProperty("CalcErrors", false,
                  "Set to true to calcuate errors when output isn't created "
                  "(default is false).");
  declareProperty("OutputCompositeMembers", false,
                  "If true and CreateOutput is true then the value of each "
                  "member of a Composite Function is also output.");
  declareProperty(new Kernel::PropertyWithValue<bool>("ConvolveMembers", false),
                  "If true and OutputCompositeMembers is true members of any "
                  "Convolution are output convolved\n"
                  "with corresponding resolution");
  declareProperty("OutputParametersOnly", false,
                  "Set to true to output only the parameters and not "
                  "workspace(s) with the calculated values\n"
                  "(default is false, ignored if CreateOutput is false and "
                  "Output is an empty string).");
  declareProperty("Surface", "", "Type of surface to output. No output if empty.");
  declareProperty("SurfaceScaling", 0.0, "Scaling factor.");
  declareProperty("SurfaceStart", 0, "Starting index.");
  declareProperty("SurfaceEnd", 0, "Ending index.");
  declareProperty(new Kernel::ArrayProperty<double>("SurfaceParams"));
}

/**
  * Copy all output workspace properties from the minimizer to Fit algorithm.
  * @param minimizer :: The minimizer to copy from.
  */
void Fit::copyMinimizerOutput(const API::IFuncMinimizer &minimizer) {
  auto &properties = minimizer.getProperties();
  for (auto prop = properties.begin(); prop != properties.end(); ++prop) {
    if ((**prop).direction() == Kernel::Direction::Output &&
        (**prop).isValid() == "") {
      Kernel::Property *property = (**prop).clone();
      declareProperty(property);
    }
  }
}

/** Executes the algorithm
*
*  @throw runtime_error Thrown if algorithm cannot execute
*/
void Fit::execConcrete() {

  std::string ties = getPropertyValue("Ties");
  if (!ties.empty()) {
    m_function->addTies(ties);
  }
  std::string contstraints = getPropertyValue("Constraints");
  if (!contstraints.empty()) {
    m_function->addConstraints(contstraints);
  }

  // prepare the function for a fit
  m_function->setUpForFit();

  API::FunctionDomain_sptr domain;
  API::FunctionValues_sptr values;
  m_domainCreator->createDomain(domain, values);

  // do something with the function which may depend on workspace
  m_domainCreator->initFunction(m_function);

  // get the minimizer
  std::string minimizerName = getPropertyValue("Minimizer");
  API::IFuncMinimizer_sptr minimizer =
      API::FuncMinimizerFactory::Instance().createMinimizer(minimizerName);

  // Try to retrieve optional properties
  int intMaxIterations = getProperty("MaxIterations");
  const size_t maxIterations = static_cast<size_t>(intMaxIterations);

  // get the cost function which must be a CostFuncFitting
  boost::shared_ptr<CostFunctions::CostFuncFitting> costFunc =
      boost::dynamic_pointer_cast<CostFunctions::CostFuncFitting>(
          API::CostFunctionFactory::Instance().create(
              getPropertyValue("CostFunction")));

  costFunc->setFittingFunction(m_function, domain, values);
  minimizer->initialize(costFunc, maxIterations);

  const int64_t nsteps = maxIterations * m_function->estimateNoProgressCalls();
  API::Progress prog(this, 0.0, 1.0, nsteps);
  m_function->setProgressReporter(&prog);
  pushParameters(*costFunc);

  // do the fitting until success or iteration limit is reached
  size_t iter = 0;
  bool success = false;
  std::string errorString;
  g_log.debug("Starting minimizer iteration\n");
  while (iter < maxIterations) {
    g_log.debug() << "Starting iteration " << iter << "\n";
    m_function->iterationStarting();
    if (!minimizer->iterate(iter)) {
      errorString = minimizer->getError();
      g_log.debug() << "Iteration stopped. Minimizer status string="
                    << errorString << "\n";

      success = errorString.empty() || errorString == "success";
      if (success) {
        errorString = "success";
      }
      break;
    }
    pushParameters(*costFunc);
    prog.report();
    m_function->iterationFinished();
    ++iter;
  }
  g_log.debug() << "Number of minimizer iterations=" << iter << "\n";

  minimizer->finalize();
  pushParameters(*costFunc);

  if (iter >= maxIterations) {
    if (!errorString.empty()) {
      errorString += '\n';
    }
    errorString += "Failed to converge after " +
                   boost::lexical_cast<std::string>(maxIterations) +
                   " iterations.";
  }

  // return the status flag
  setPropertyValue("OutputStatus", errorString);

  // degrees of freedom
  size_t dof = domain->size() - costFunc->nParams();
  if (dof == 0)
    dof = 1;
  double rawcostfuncval = minimizer->costFunctionVal();
  double finalCostFuncVal = rawcostfuncval / double(dof);

  setProperty("OutputChi2overDoF", finalCostFuncVal);

  // fit ended, creating output

  // get the workspace
  API::Workspace_const_sptr ws = getProperty("InputWorkspace");

  bool doCreateOutput = getProperty("CreateOutput");
  std::string baseName = getPropertyValue("Output");
  if (!baseName.empty()) {
    doCreateOutput = true;
  }
  bool doCalcErrors = getProperty("CalcErrors");
  if (doCreateOutput) {
    doCalcErrors = true;
  }
  if (costFunc->nParams() == 0) {
    doCalcErrors = false;
  }

  GSLMatrix covar;
  if (doCalcErrors) {
    // Calculate the covariance matrix and the errors.
    costFunc->calCovarianceMatrix(covar);
    costFunc->calFittingErrors(covar, rawcostfuncval);
  }

  if (doCreateOutput) {
    copyMinimizerOutput(*minimizer);

    if (baseName.empty()) {
      baseName = ws->name();
      if (baseName.empty()) {
        baseName = "Output";
      }
    }
    baseName += "_";

    declareProperty(
        new API::WorkspaceProperty<API::ITableWorkspace>(
            "OutputNormalisedCovarianceMatrix", "", Kernel::Direction::Output),
        "The name of the TableWorkspace in which to store the final covariance "
        "matrix");
    setPropertyValue("OutputNormalisedCovarianceMatrix",
                     baseName + "NormalisedCovarianceMatrix");

    Mantid::API::ITableWorkspace_sptr covariance =
        Mantid::API::WorkspaceFactory::Instance().createTable("TableWorkspace");
    covariance->addColumn("str", "Name");
    // set plot type to Label = 6
    covariance->getColumn(covariance->columnCount() - 1)->setPlotType(6);
    // std::vector<std::string> paramThatAreFitted; // used for populating 1st
    // "name" column
    for (size_t i = 0; i < m_function->nParams(); i++) {
      if (m_function->isActive(i)) {
        covariance->addColumn("double", m_function->parameterName(i));
        // paramThatAreFitted.push_back(m_function->parameterName(i));
      }
    }

    size_t np = m_function->nParams();
    size_t ia = 0;
    for (size_t i = 0; i < np; i++) {
      if (m_function->isFixed(i))
        continue;
      Mantid::API::TableRow row = covariance->appendRow();
      row << m_function->parameterName(i);
      size_t ja = 0;
      for (size_t j = 0; j < np; j++) {
        if (m_function->isFixed(j))
          continue;
        if (j == i)
          row << 100.0;
        else {
          if (!covar.gsl()) {
            throw std::runtime_error(
                "There was an error while allocating the (GSL) covariance "
                "matrix "
                "which is needed to produce fitting error results.");
          }
          row << 100.0 * covar.get(ia, ja) /
                     sqrt(covar.get(ia, ia) * covar.get(ja, ja));
        }
        ++ja;
      }
      ++ia;
    }

    setProperty("OutputNormalisedCovarianceMatrix", covariance);

    // create output parameter table workspace to store final fit parameters
    // including error estimates if derivative of fitting function defined

    declareProperty(new API::WorkspaceProperty<API::ITableWorkspace>(
                        "OutputParameters", "", Kernel::Direction::Output),
                    "The name of the TableWorkspace in which to store the "
                    "final fit parameters");

    setPropertyValue("OutputParameters", baseName + "Parameters");

    Mantid::API::ITableWorkspace_sptr result =
        Mantid::API::WorkspaceFactory::Instance().createTable("TableWorkspace");
    result->addColumn("str", "Name");
    // set plot type to Label = 6
    result->getColumn(result->columnCount() - 1)->setPlotType(6);
    result->addColumn("double", "Value");
    result->addColumn("double", "Error");
    // yErr = 5
    result->getColumn(result->columnCount() - 1)->setPlotType(5);

    for (size_t i = 0; i < m_function->nParams(); i++) {
      Mantid::API::TableRow row = result->appendRow();
      row << m_function->parameterName(i) << m_function->getParameter(i)
          << m_function->getError(i);
    }
    // Add chi-squared value at the end of parameter table
    Mantid::API::TableRow row = result->appendRow();
#if 1
    std::string costfuncname = getPropertyValue("CostFunction");
    if (costfuncname == "Rwp")
      row << "Cost function value" << rawcostfuncval;
    else
      row << "Cost function value" << finalCostFuncVal;
    setProperty("OutputParameters", result);
#else
    row << "Cost function value" << finalCostFuncVal;
    Mantid::API::TableRow row2 = result->appendRow();
    std::string name(getPropertyValue("CostFunction"));
    name += " value";
    row2 << name << rawcostfuncval;
#endif

    setProperty("OutputParameters", result);

    bool outputParametersOnly = getProperty("OutputParametersOnly");

    if (!outputParametersOnly) {
      const bool unrollComposites = getProperty("OutputCompositeMembers");
      bool convolveMembers = existsProperty("ConvolveMembers");
      if (convolveMembers) {
        convolveMembers = getProperty("ConvolveMembers");
      }
      m_domainCreator->separateCompositeMembersInOutput(unrollComposites,
                                                        convolveMembers);
      m_domainCreator->createOutputWorkspace(baseName, m_function, domain,
                                             values);
    }
  }

  outputSurface();

  progress(1.0);
}

void Fit::pushParameters(const CostFunctions::CostFuncFitting& fun) {
  GSLVector params;
  fun.getParameters(params);
  m_points.push_back(params);
}

void Fit::outputSurface() {
  std::string surface = getPropertyValue("Surface");
  if (surface.empty()) return;

  size_t n = 100;
  auto matrix = 
    Mantid::API::WorkspaceFactory::Instance().create("Workspace2D", n, n, n);
  API::Axis *const verticalAxis = new API::NumericAxis(n);
  //verticalAxis->unit() =
  //    Kernel::UnitFactory::Instance().create("Label");
  matrix->replaceAxis(1, verticalAxis);

  std::string outputBaseName = getPropertyValue("Output");
  if (outputBaseName.empty()) {
    outputBaseName = "out";
  }

  try {
    int ii0 = getProperty("SurfaceStart");
    int ii1 = getProperty("SurfaceEnd");
    size_t i0 = static_cast<size_t>(ii0);
    if (i0 >= m_points.size() - 1) {
      i0 = 0;
    }
    size_t i1 = static_cast<size_t>(ii1);
    if (i1 == 0 || i1 >= m_points.size() - 1) {
      i1 = m_points.size() - 1;
    }
    const GSLVector& p0 = m_points[i0];
    GSLVector p1 = m_points[i1];

    std::vector<double> finalParams = getProperty("SurfaceParams");
    if (!finalParams.empty()) {
      if (p1.size() != finalParams.size()) {
        std::cerr << "Wrong vector size of SurfaceParams" << std::endl;
      }
      for(size_t i = 0; i < p1.size(); ++i) {
        p1[i] = finalParams[i];
      }
    }

    auto funCopy = m_function->clone();

    std::cerr << m_points.size() << " iterations" << std::endl;

    API::FunctionDomain_sptr domain;
    API::FunctionValues_sptr values;
    m_domainCreator->createDomain(domain, values);
    m_domainCreator->initFunction(funCopy);
    boost::shared_ptr<CostFunctions::CostFuncFitting> costFunc =
        boost::dynamic_pointer_cast<CostFunctions::CostFuncFitting>(
            API::CostFunctionFactory::Instance().create(
                getPropertyValue("CostFunction")));

    costFunc->setFittingFunction(funCopy, domain, values);
    FuncMinimisers::LocalSearchMinimizer minimizer;
    minimizer.initialize(costFunc);
    auto val1 = costFunc->val();
    costFunc->setParameters(p0);
    auto val0 = costFunc->val();
    costFunc->setParameters(p1);

    auto y = p1;
    y -= p0;
    double yWidth = y.norm() * 2;
    y.normalize();
    double dy = yWidth / static_cast<double>(n);

    auto x = costFunc->getDeriv();
    {
      x *= -1;
      auto dotProd = x.dot(y);
      auto pd = y;
      pd *= dotProd;
      x -= pd;
      x.normalize();
    }
    std::cerr << "p0 " << p0 << std::endl;
    std::cerr << "p1 " << p1 << std::endl;

    std::cerr << "x " << x << std::endl;
    std::cerr << "y " << y << std::endl;
    std::cerr << x.dot(y) << std::endl;

    double dx = yWidth / 2;
    {
      const double shift = 1e-6;
      auto dp = p1;
      dp *= shift;
      for(size_t i = 0; i < dp.size(); ++i) {
        if (dp[i] == 0.0) {
          dp[i] = shift;
        }
      }
      std::cerr << "dp " << dp << std::endl;
      dx = dp.norm();
    }

    declareProperty(
        new API::WorkspaceProperty<API::ITableWorkspace>(
            "SurfacePath", "", Kernel::Direction::Output),
        "The name of the TableWorkspace ");
    setPropertyValue("SurfacePath", outputBaseName + "_Path");
    Mantid::API::ITableWorkspace_sptr path =
        Mantid::API::WorkspaceFactory::Instance().createTable("TableWorkspace");
    path->addColumn("double", "X")->setPlotType(1);
    path->addColumn("double", "Y")->setPlotType(2);

    for (size_t i = i0; i <= i1; i++) {
      auto p = m_points[i];
      p -= p0;
      auto px = p.dot(x);
      if (fabs(px) > dx) {
        dx = fabs(px);
      }
      Mantid::API::TableRow row = path->appendRow();
      row << px << p.dot(y) - yWidth / 2;
    }
    dx *= 2;
    setProperty("SurfacePath", path);

    double scaling = getProperty("SurfaceScaling");
    if (scaling != 0.0) {
      dx *= scaling;
    } else {
      double oldVal = val1;
      for (size_t j = 0; j < 100; ++j) {
        auto px = x;
        px *= dx;
        px += p1;
        costFunc->setParameters(px);
        auto val = costFunc->val();
        auto ratio = fabs(val / val0);
        if (j == 0 && ratio >= 1.0) {
          break;
        }
        if (j > 0 && ratio < 1.1 && ratio > 0.9) {
          if (val != oldVal) {
            std::cerr << "Search stopped after " << j << " iterations " << val0
                      << ' ' << val << std::endl;
            break;
          }
        }
        if (val < val0) {
          dx *= 2;
        } else {
          dx *= 0.75;
        }
        oldVal = val;
        if (j == 99) {
          std::cerr << "Search didn't converge" << std::endl;
        }
      }
    }
    dx /= n;

    std::cerr << "dx=" << dx << std::endl;
    std::cerr << "dy=" << dy << std::endl;
    for (size_t i = 0; i < n; ++i) {
      double yi = i * dy;
      auto py = y;
      py *= yi;
      verticalAxis->setValue(i, yi - yWidth / 2);
      for (size_t j = 0; j < n; ++j) {
        double xj = j * dx - dx * n / 2;
        auto px = x;
        px *= xj;
        px += py;
        px += p0;
        costFunc->setParameters(px);
        auto val = costFunc->val();
        matrix->dataX(i)[j] = xj;
        matrix->dataY(i)[j] = val;
      }
    }

  } catch (std::runtime_error &e) {
    g_log.warning() << "Cannot create the surface, error happened:" << std::endl
                    << e.what() << std::endl;
  }

  declareProperty(
      new API::WorkspaceProperty<API::MatrixWorkspace>(
          "SurfaceWorkspace", "", Kernel::Direction::Output),
      "Name of the output Workspace holding resulting simulated spectrum");
  setPropertyValue("SurfaceWorkspace", outputBaseName + "_Surface");
  setProperty("SurfaceWorkspace", matrix);
}

} // namespace Algorithms
} // namespace CurveFitting
} // namespace Mantid
