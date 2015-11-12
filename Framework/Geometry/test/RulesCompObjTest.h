#ifndef MANTID_RULESCOMPOBJTEST__
#define MANTID_RULESCOMPOBJTEST__
#include <cxxtest/TestSuite.h>
#include <cmath>
#include <vector>
#include "MantidKernel/Logger.h"
#include "MantidKernel/System.h"
#include <cfloat>
#include "MantidKernel/V3D.h"
#include "MantidGeometry/Surfaces/Quadratic.h"
#include "MantidGeometry/Objects/Object.h"
#include "MantidGeometry/Objects/Rules.h"
#include "MantidGeometry/Surfaces/Plane.h"
#include "MantidGeometry/Surfaces/Sphere.h"
#include "MantidGeometry/Surfaces/Cylinder.h"
#include "MantidGeometry/Surfaces/Cone.h"

#include "boost/make_shared.hpp"

using namespace Mantid;
using namespace Geometry;
using Mantid::Kernel::V3D;

class RulesCompObjTest : public CxxTest::TestSuite {
public:
  void testConstructor() {
    CompObj A;
    TS_ASSERT_EQUALS(A.leaf(0), (Rule *)0);
    TS_ASSERT_EQUALS(A.leaf(1), (Rule *)0);
    TS_ASSERT_EQUALS(A.display(), "#0");
    TS_ASSERT_EQUALS(A.getObjN(), 0);
    TS_ASSERT_EQUALS(A.getObj(), (Object *)0);
    TS_ASSERT_EQUALS(A.isComplementary(), 1);
  }

  void testSetObject() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    TS_ASSERT_EQUALS(A.display(), "#10");
    TS_ASSERT_EQUALS(A.getObjN(), 10);
    TS_ASSERT_EQUALS(A.getObj(), &cpCylinder);
  }

  void testCompObjConstructor() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj B(A);
    TS_ASSERT_EQUALS(B.display(), "#10");
    TS_ASSERT_EQUALS(B.getObjN(), 10);
    TS_ASSERT_EQUALS(B.getObj(), &cpCylinder);
  }

  void testClone() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj *B;
    B = A.clone();
    TS_ASSERT_EQUALS(B->display(), "#10");
    TS_ASSERT_EQUALS(B->getObjN(), 10);
    TS_ASSERT_EQUALS(B->getObj(), &cpCylinder);
    delete B;
  }

  void testAssignment() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj B;
    B = A;
    TS_ASSERT_EQUALS(B.display(), "#10");
    TS_ASSERT_EQUALS(B.getObjN(), 10);
    TS_ASSERT_EQUALS(B.getObj(), &cpCylinder);
  }

  void testSetLeaves() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj B;
    B.setLeaves(&A, (Rule *)0);
    TS_ASSERT_EQUALS(B.display(), "#10");
    TS_ASSERT_EQUALS(B.getObjN(), 10);
    TS_ASSERT_EQUALS(B.getObj(), &cpCylinder);
  }

  void testSetLeaf() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj B;
    B.setLeaf(&A, 0);
    TS_ASSERT_EQUALS(B.display(), "#10");
    TS_ASSERT_EQUALS(B.getObjN(), 10);
    TS_ASSERT_EQUALS(B.getObj(), &cpCylinder);
  }

  void testFindLeaf() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj B;
    TS_ASSERT_EQUALS(A.findLeaf(&A), 0);
    TS_ASSERT_EQUALS(A.findLeaf(&B), -1);
    B = A;
    TS_ASSERT_EQUALS(A.findLeaf(&B), -1);
  }

  void testFindKey() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    CompObj B;
    TS_ASSERT_EQUALS(A.findKey(10), (Rule *)0); // Always returns 0
    TS_ASSERT_EQUALS(A.findKey(11), (Rule *)0);
  }

  void testIsValid() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    TS_ASSERT_EQUALS(
        A.isValid(V3D(0.0, 0.0, 0.0)),
        false); // center is inside the cylinder so it will return complement
    TS_ASSERT_EQUALS(A.isValid(V3D(1.3, 0.0, 0.0)),
                     true); // outside cap cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(1.2, 0.0, 0.0)),
                     false); // on the cap end of cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(1.1, 0.0, 0.0)),
                     false); // inside the cap end of cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(-3.3, 0.0, 0.0)),
                     true); // outside other end of cap cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(-3.2, 0.0, 0.0)),
                     false); // on end of cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(-3.1, 0.0, 0.0)),
                     false); // inside the cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(0.0, 3.1, 0.0)), true);  // outside cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(0.0, 3.0, 0.0)), false); // on the cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(0.0, 2.9, 0.0)), false); // inside cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(0.0, 0.0, 3.1)), true);  // outside cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(0.0, 0.0, 3.0)), false); // on the cylinder
    TS_ASSERT_EQUALS(A.isValid(V3D(0.0, 0.0, 2.9)), false); // inside cylinder
  }

  void testIsValidMap() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);

    std::map<int, int> input;
    input[31] = 1;
    input[32] = 1;
    input[33] = 1;
    TS_ASSERT_EQUALS(A.isValid(input), true);
    input[31] = 0;
    TS_ASSERT_EQUALS(A.isValid(input), true);
    input[32] = 0;
    TS_ASSERT_EQUALS(A.isValid(input), false);
    input[33] = 0;
    TS_ASSERT_EQUALS(A.isValid(input), true);
    input[32] = 1;
    TS_ASSERT_EQUALS(A.isValid(input), true);
    input[33] = 1;
    TS_ASSERT_EQUALS(A.isValid(input), true);
  }

  void testSimplyfy() {
    Object cpCylinder = createCappedCylinder();
    CompObj A;
    A.setObj(&cpCylinder);
    A.setObjN(10);
    TS_ASSERT_EQUALS(A.simplify(),
                     0); // Always return 0 bcos a end node cannot be simplified
  }

private:
  Object createCappedCylinder() {
    std::string C31 = "cx 3.0"; // cylinder x-axis radius 3
    std::string C32 = "px 1.2";
    std::string C33 = "px -3.2";

    // First create some surfaces
    std::map<int, boost::shared_ptr<Surface>> CylSurMap;
    CylSurMap[31] = boost::make_shared<Cylinder>();
    CylSurMap[32] = boost::make_shared<Plane>();
    CylSurMap[33] = boost::make_shared<Plane>();

    CylSurMap[31]->setSurface(C31);
    CylSurMap[32]->setSurface(C32);
    CylSurMap[33]->setSurface(C33);
    CylSurMap[31]->setName(31);
    CylSurMap[32]->setName(32);
    CylSurMap[33]->setName(33);

    // Capped cylinder (id 21)
    // using surface ids: 31 (cylinder) 32 (plane (top) ) and 33 (plane (base))
    std::string ObjCapCylinder = "-31 -32 33";

    Object retVal;
    retVal.setObject(21, ObjCapCylinder);
    retVal.populate(CylSurMap);

    return retVal;
  }
};
//---------------------------------End of
// CompObj----------------------------------------

#endif