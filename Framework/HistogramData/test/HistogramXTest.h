#ifndef MANTID_HISTOGRAMDATA_HISTOGRAMXTEST_H_
#define MANTID_HISTOGRAMDATA_HISTOGRAMXTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidHistogramData/HistogramX.h"

//using Mantid::HistogramData::Points;
//using Mantid::HistogramData::BinEdges;
using Mantid::HistogramData::HistogramX;

class HistogramXTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static HistogramXTest *createSuite() { return new HistogramXTest(); }
  static void destroySuite(HistogramXTest *suite) { delete suite; }

  void test_fake() {}
  /*
  void xtest_construct_from_Points() {
    Points points{0.1, 0.2, 0.4};
    const HistogramX histX(points);
    TS_ASSERT_EQUALS(histX.size(), 3);
  }

  void xtest_construct_from_BinEdges() {
    BinEdges binEdges{0.1, 0.2, 0.4};
    const HistogramX histX(binEdges);
    TS_ASSERT_EQUALS(histX.size(), 3);
  }

  void xtest_construct_from_invalid_BinEdges() {
    BinEdges binEdges(1);
    TS_ASSERT_THROWS(HistogramX{binEdges}, std::logic_error);
  }

  void xtest_assignment() {
    const HistogramX src(Points(1));
    HistogramX dest(Points(1));
    TS_ASSERT_THROWS_NOTHING(dest = src);
    const auto &constDest = static_cast<const HistogramX>(dest);
    TS_ASSERT_EQUALS(&constDest[0], &src[0]);
  }

  void xtest_assignment_mutating() {
    const HistogramX src(Points(1));
    HistogramX dest(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(dest = src);
    const auto &constDest = static_cast<const HistogramX>(dest);
    TS_ASSERT_EQUALS(&constDest[0], &src[0]);
  }

  void xtest_assignment_failures() {
    HistogramX src1(Points(1));
    HistogramX dest1(Points(2));
    TS_ASSERT_THROWS(dest1 = src1, std::logic_error);
    HistogramX src2(Points(1));
    HistogramX dest2(BinEdges(3));
    TS_ASSERT_THROWS(dest2 = src2, std::logic_error);
    HistogramX src3(BinEdges(2));
    HistogramX dest3(Points(2));
    TS_ASSERT_THROWS(dest3 = src3, std::logic_error);
    HistogramX src4(BinEdges(2));
    HistogramX dest4(BinEdges(3));
    TS_ASSERT_THROWS(dest4 = src4, std::logic_error);
  }

  void xtest_points_from_edges() {
    BinEdges binEdges{0.1, 0.2, 0.4};
    const HistogramX histX(binEdges);
    TS_ASSERT_EQUALS(histX.size(), 3);
    const auto points = histX.points();
    TS_ASSERT_DIFFERS(&points[0], &histX[0]);
    TS_ASSERT_EQUALS(points.size(), 2);
    TS_ASSERT_DELTA(points[0], 0.15, 1e-14);
    TS_ASSERT_DELTA(points[1], 0.3, 1e-14);
  }

  void xtest_points_from_points() {
    const HistogramX histX(Points{0.1, 0.2, 0.4});
    const auto points = histX.points();
    TS_ASSERT_EQUALS(&points[0], &histX[0]);
  }

  void xtest_setPoints_from_vector() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(x1.setPoints(std::vector<double>{0.1, 0.2}));
    TS_ASSERT_EQUALS(x1.size(), 2);
    TS_ASSERT_EQUALS(x1[0], 0.1);
    TS_ASSERT_EQUALS(x1[1], 0.2);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setPoints(std::vector<double>{0.1}));
    TS_ASSERT_EQUALS(x2.size(), 1);
    TS_ASSERT_EQUALS(x2[0], 0.1);
  }

  void xtest_setPoints_from_Points() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(x1.setPoints(Points{0.1, 0.2}));
    TS_ASSERT_EQUALS(x1.size(), 2);
    TS_ASSERT_EQUALS(x1[0], 0.1);
    TS_ASSERT_EQUALS(x1[1], 0.2);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setPoints(Points{0.1}));
    TS_ASSERT_EQUALS(x2.size(), 1);
    TS_ASSERT_EQUALS(x2[0], 0.1);
  }

  void xtest_setPoints_from_BinEdges() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(x1.setPoints(BinEdges{0.1, 0.2, 0.4}));
    TS_ASSERT_EQUALS(x1.size(), 2);
    TS_ASSERT_DELTA(x1[0], 0.15, 1e-14);
    TS_ASSERT_DELTA(x1[1], 0.3, 1e-14);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setPoints(BinEdges{0.1, 0.2}));
    TS_ASSERT_EQUALS(x2.size(), 1);
    TS_ASSERT_DELTA(x2[0], 0.15, 1e-14);
  }

  void xtest_setPoints_degenerate() {
    HistogramX x1(Points(0));
    TS_ASSERT_THROWS_NOTHING(x1.setPoints(std::vector<double>(0)));
    TS_ASSERT_EQUALS(x1.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x1.setPoints(Points(0)));
    TS_ASSERT_EQUALS(x1.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x1.setPoints(BinEdges(0)));
    TS_ASSERT_EQUALS(x1.size(), 0);
    HistogramX x2(BinEdges(0));
    TS_ASSERT_THROWS_NOTHING(x2.setPoints(std::vector<double>(0)));
    TS_ASSERT_EQUALS(x2.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x2.setPoints(Points(0)));
    TS_ASSERT_EQUALS(x2.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x2.setPoints(BinEdges(0)));
    TS_ASSERT_EQUALS(x2.size(), 0);
  }

  void xtest_setPoints_size_mismatch() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS(x1.setPoints(std::vector<double>(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(std::vector<double>(3)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(Points(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(Points(3)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(BinEdges(2)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(BinEdges(4)), std::logic_error);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS(x2.setPoints(std::vector<double>(0)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(std::vector<double>(2)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(Points(0)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(Points(2)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(BinEdges(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(BinEdges(3)), std::logic_error);
  }

  void xtest_setPoints_size_mismatch_degenerate() {
    HistogramX x1(Points(0));
    TS_ASSERT_THROWS(x1.setPoints(std::vector<double>(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(Points(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setPoints(BinEdges(1)), std::logic_error);
    HistogramX x2(BinEdges(0));
    TS_ASSERT_THROWS(x2.setPoints(std::vector<double>(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(Points(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setPoints(BinEdges(1)), std::logic_error);
  }

  void xtest_edges_from_edges() {
    const HistogramX histX(BinEdges{0.1, 0.2, 0.4});
    const auto edges = histX.binEdges();
    TS_ASSERT_EQUALS(&edges[0], &histX[0]);
    TS_ASSERT_EQUALS(edges.size(), 3);
  }

  void xtest_edges_from_points() {
    const HistogramX histX(Points{0.1, 0.2, 0.4});
    const auto edges = histX.binEdges();
    TS_ASSERT_DIFFERS(&edges[0], &histX[0]);
    TS_ASSERT_EQUALS(edges.size(), 4);
  }

  void xtest_setBinEdges() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(std::vector<double>(3)));
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(Points(2)));
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(BinEdges(3)));
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(std::vector<double>(2)));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(Points(1)));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(BinEdges(2)));
  }

  void xtest_setBinEdges_from_vector() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(
        x1.setBinEdges(std::vector<double>{0.1, 0.2, 0.4}));
    TS_ASSERT_EQUALS(x1.size(), 3);
    TS_ASSERT_EQUALS(x1[0], 0.1);
    TS_ASSERT_EQUALS(x1[1], 0.2);
    TS_ASSERT_EQUALS(x1[2], 0.4);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(std::vector<double>{0.1, 0.2}));
    TS_ASSERT_EQUALS(x2.size(), 2);
    TS_ASSERT_EQUALS(x2[0], 0.1);
    TS_ASSERT_EQUALS(x2[1], 0.2);
  }

  void xtest_setBinEdges_from_Points() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(Points{0.1, 0.3}));
    TS_ASSERT_EQUALS(x1.size(), 3);
    TS_ASSERT_DELTA(x1[0], 0.0, 1e-14);
    TS_ASSERT_DELTA(x1[1], 0.2, 1e-14);
    TS_ASSERT_DELTA(x1[2], 0.4, 1e-14);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(Points{1.0}));
    TS_ASSERT_EQUALS(x2.size(), 2);
    TS_ASSERT_DELTA(x2[0], 0.5, 1e-14);
    TS_ASSERT_DELTA(x2[1], 1.5, 1e-14);
  }

  void xtest_setBinEdges_from_BinEdges() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(BinEdges{0.1, 0.2, 0.4}));
    TS_ASSERT_EQUALS(x1.size(), 3);
    TS_ASSERT_EQUALS(x1[0], 0.1);
    TS_ASSERT_EQUALS(x1[1], 0.2);
    TS_ASSERT_EQUALS(x1[2], 0.4);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(BinEdges{0.1, 0.2}));
    TS_ASSERT_EQUALS(x2.size(), 2);
    TS_ASSERT_EQUALS(x2[0], 0.1);
    TS_ASSERT_EQUALS(x2[1], 0.2);
  }

  void xtest_setBinEdges_degenerate() {
    HistogramX x1(Points(0));
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(std::vector<double>(0)));
    TS_ASSERT_EQUALS(x1.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(Points(0)));
    TS_ASSERT_EQUALS(x1.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x1.setBinEdges(BinEdges(0)));
    TS_ASSERT_EQUALS(x1.size(), 0);
    HistogramX x2(BinEdges(0));
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(std::vector<double>(0)));
    TS_ASSERT_EQUALS(x2.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(Points(0)));
    TS_ASSERT_EQUALS(x2.size(), 0);
    TS_ASSERT_THROWS_NOTHING(x2.setBinEdges(BinEdges(0)));
    TS_ASSERT_EQUALS(x2.size(), 0);
  }

  void xtest_setBinEdges_size_mismatch() {
    HistogramX x1(Points(2));
    TS_ASSERT_THROWS(x1.setBinEdges(std::vector<double>(2)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(std::vector<double>(4)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(Points(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(Points(3)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(BinEdges(2)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(BinEdges(4)), std::logic_error);
    HistogramX x2(BinEdges(2));
    TS_ASSERT_THROWS(x2.setBinEdges(std::vector<double>(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(std::vector<double>(3)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(Points(0)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(Points(2)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(BinEdges(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(BinEdges(3)), std::logic_error);
  }

  void xtest_setBinEdges_size_mismatch_degenerate() {
    HistogramX x1(Points(0));
    TS_ASSERT_THROWS(x1.setBinEdges(std::vector<double>(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(Points(1)), std::logic_error);
    TS_ASSERT_THROWS(x1.setBinEdges(BinEdges(1)), std::logic_error);
    HistogramX x2(BinEdges(0));
    TS_ASSERT_THROWS(x2.setBinEdges(std::vector<double>(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(Points(1)), std::logic_error);
    TS_ASSERT_THROWS(x2.setBinEdges(BinEdges(1)), std::logic_error);
  }
  */
};

#endif /* MANTID_HISTOGRAMDATA_HISTOGRAMXTEST_H_ */
