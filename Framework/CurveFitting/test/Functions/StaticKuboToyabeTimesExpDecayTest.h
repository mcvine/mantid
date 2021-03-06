#ifndef MANTID_CURVEFITTING_STATICKUBOTOYABETIMESEXPDECAYTEST_H_
#define MANTID_CURVEFITTING_STATICKUBOTOYABETIMESEXPDECAYTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidCurveFitting/Functions/StaticKuboToyabeTimesExpDecay.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidCurveFitting/Algorithms/Fit.h"
#include "MantidDataObjects/Workspace2D.h"

using Mantid::CurveFitting::Functions::StaticKuboToyabeTimesExpDecay;

using namespace Mantid::Kernel;
using namespace Mantid::API;
using namespace Mantid::CurveFitting;
using namespace Mantid::CurveFitting::Functions;
using namespace Mantid::CurveFitting::Algorithms;
using namespace Mantid::DataObjects;

class StaticKuboToyabeTimesExpDecayTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static StaticKuboToyabeTimesExpDecayTest *createSuite() {
    return new StaticKuboToyabeTimesExpDecayTest();
  }
  static void destroySuite(StaticKuboToyabeTimesExpDecayTest *suite) {
    delete suite;
  }

  void getMockData(Mantid::MantidVec &y, Mantid::MantidVec &e) {
    // A = 0.24, Delta = 0.16, Lambda = 0.1
    y[0] = 0.24;
    y[1] = 0.211661;
    y[2] = 0.177213;
    y[3] = 0.140561;
    y[4] = 0.10522;
    y[5] = 0.0738913;
    y[6] = 0.0482474;
    y[7] = 0.0289314;
    y[8] = 0.015716;
    y[9] = 0.00776156;
    y[10] = 0.00390022;
    y[11] = 0.00288954;
    y[12] = 0.00360064;
    y[13] = 0.00512831;
    y[14] = 0.00682993;

    for (int i = 0; i < 15; i++)
      e[i] = 1.0;
  }

  StaticKuboToyabeTimesExpDecayTest() : fn() {}

  void test_Initialize() { TS_ASSERT_THROWS_NOTHING(fn.initialize()); }

  void test_Name() {
    TS_ASSERT_EQUALS(fn.name(), "StaticKuboToyabeTimesExpDecay");
  }

  void test_Params() {
    TS_ASSERT_DELTA(fn.getParameter("A"), 0.2, 0.0001);
    TS_ASSERT_DELTA(fn.getParameter("Delta"), 0.2, 0.0001);
    TS_ASSERT_DELTA(fn.getParameter("Lambda"), 0.2, 0.0001);
  }

  void test_Category() {
    const std::vector<std::string> categories = fn.categories();
    TS_ASSERT(categories.size() == 1);
    TS_ASSERT(categories[0] == "Muon");
  }

  void test_AgainstMockData() {
    Algorithms::Fit alg2;
    TS_ASSERT_THROWS_NOTHING(alg2.initialize());
    TS_ASSERT(alg2.isInitialized());

    // create mock data to test against
    std::string wsName = "SKTTimesExpDecayMockData";
    Workspace_sptr ws =
        WorkspaceFactory::Instance().create("Workspace2D", 1, 15, 15);
    Workspace2D_sptr ws2D = boost::dynamic_pointer_cast<Workspace2D>(ws);

    for (int i = 0; i < 15; i++)
      ws2D->dataX(0)[i] = i;

    getMockData(ws2D->dataY(0), ws2D->dataE(0));

    // put this workspace in the data service
    TS_ASSERT_THROWS_NOTHING(
        AnalysisDataService::Instance().addOrReplace(wsName, ws2D));

    alg2.setPropertyValue("Function", fn.asString());

    // Set which spectrum to fit against and initial starting values
    alg2.setPropertyValue("InputWorkspace", wsName);
    alg2.setPropertyValue("WorkspaceIndex", "0");
    alg2.setPropertyValue("StartX", "0");
    alg2.setPropertyValue("EndX", "14");

    TS_ASSERT_THROWS_NOTHING(TS_ASSERT(alg2.execute()))

    TS_ASSERT(alg2.isExecuted());

    double dummy = alg2.getProperty("OutputChi2overDoF");
    TS_ASSERT_DELTA(dummy, 0.0001, 0.0001);

    IFunction_sptr out = alg2.getProperty("Function");
    TS_ASSERT_DELTA(out->getParameter("A"), 0.24, 0.0001);
    TS_ASSERT_DELTA(out->getParameter("Delta"), 0.16, 0.001);
    TS_ASSERT_DELTA(out->getParameter("Lambda"), 0.1, 0.001);

    AnalysisDataService::Instance().remove(wsName);
  }

  StaticKuboToyabeTimesExpDecay fn;
};

#endif /* MANTID_CURVEFITTING_STATICKUBOTOYABETIMESEXPDECAYTEST_H_ */