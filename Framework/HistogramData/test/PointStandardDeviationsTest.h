#ifndef MANTID_HISTOGRAMDATA_POINTSTANDARDDEVIATIONSTEST_H_
#define MANTID_HISTOGRAMDATA_POINTSTANDARDDEVIATIONSTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidHistogramData/PointStandardDeviations.h"

using Mantid::HistogramData::PointStandardDeviations;

class PointStandardDeviationsTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static PointStandardDeviationsTest *createSuite() {
    return new PointStandardDeviationsTest();
  }
  static void destroySuite(PointStandardDeviationsTest *suite) { delete suite; }

  void test_construct_default() {
    const PointStandardDeviations points{};
    TS_ASSERT(!points);
  }
};

#endif /* MANTID_HISTOGRAMDATA_POINTSTANDARDDEVIATIONSTEST_H_ */
