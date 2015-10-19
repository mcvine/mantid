#ifndef MANTID_GEOMETRY_REFLECTIONGENERATOR_H_
#define MANTID_GEOMETRY_REFLECTIONGENERATOR_H_

#include "MantidGeometry/DllConfig.h"
#include "MantidGeometry/Crystal/CrystalStructure.h"
#include "MantidGeometry/Crystal/StructureFactorCalculator.h"
#include "MantidGeometry/Crystal/HKLFilter.h"

namespace Mantid {
namespace Geometry {

enum struct ReflectionConditionFilter {
  None,
  Centering,
  SpaceGroup,
  StructureFactor
};

/** ReflectionGenerator

  ReflectionGenerator is a class that provides the means to perform some
  common tasks involving generation of reflections. While the combination
  of HKLGenerator and HKLFilter is very flexible, very often a limited set
  of operations has to be performed repeatedly, involving the crystal
  structure.

  ReflectionGenerator is constructed from a CrystalStructure object, which
  is then stored internally. Additionally, a default filter for reflection
  conditions can be set, which is applied for HKL-generation in addition
  to a DRangeFilter. For more flexibility, a method is provided that accepts
  an HKLFilter as additional argument, this filter is then added to the
  DRangeFilter.

  This way it's very simple to obtain for example a list of unique reflections
  for a given crystal structure:

    CrystalStructure structure("5.43 5.43 5.43",
                               "F d -3 m", "Si 0 0 0 1.0 0.05");
    ReflectionGenerator generator(structure);

    // Get all unique HKLs between 0.5 and 5.0 Angstrom
    std::vector<V3D> hkls = generator.getUniqueHKLs(0.5, 5.0);

  Additionally there are methods to obtain structure factors and d-values
  for a given list of HKLs.

  The generated reflection lists can be filtered by several criteria to remove
  reflections that are not allowed. In the default case, the reflection
  conditions of the space group are used. Alternatively, the reflections can be
  filtered according to the centering of the lattice or, with some more
  computational effort, by their structure factors.

  The default filter method can be supplied to the constructor in the form of
  the enum ReflectionConditionFilter. Furthermore, it's possible to provide
  self-defined filters in overloaded versions of the HKL generation methods.
  These can also be generated by using a small helper method that creates
  the filters and populates them with values associated to the stored crystal
  structure:


    ReflectionGenerator generator(
                            CrystalStructure("5.43 5.43 5.43",
                                    "F d -3 m", "Si 0 0 0 1.0 0.05"));
    auto filter = generator.getReflectionConditionFilter(
                                    ReflectionConditionFilter::StructureFactor);

  The filter can then be used in connection with HKLGenerator or the HKL-
  generation methods of ReflectionGenerator.

  An example where structure factors are required can be found in a very
  common crystal structure, the structure of Silicon. Silicon crystallizes
  in the space group Fd-3m (No. 227) and the asymmetric unit consists of
  one Si-atom at the position (1/8, 1/8, 1/8) (for origin choice 2, with the
  inversion at the origin). Looking up this space group in the International
  Tables for Crystallography A reveals that placing a scatterer at this
  position introduces a new reflection condition for general reflections hkl:

         h = 2n + 1 (reflections with odd h are allowed)
      or h + k + l = 4n

  This means that for example the reflection family {2 2 2} is not allowed,
  even though the F-centering would allow it. Using the structure factor
  calculation filter solves this problem.

      @author Michael Wedel, ESS
      @date 30/09/2015

  Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
  National Laboratory & European Spallation Source

  This file is part of Mantid.

  Mantid is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  Mantid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  File change history is stored at: <https://github.com/mantidproject/mantid>
  Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class MANTID_GEOMETRY_DLL ReflectionGenerator {
public:
  ReflectionGenerator(const CrystalStructure &crystalStructure,
                      ReflectionConditionFilter defaultFilter =
                          ReflectionConditionFilter::SpaceGroup);
  ~ReflectionGenerator() {}

  const CrystalStructure &getCrystalStructure() const;

  HKLFilter_const_sptr getDRangeFilter(double dMin, double dMax) const;
  HKLFilter_const_sptr
  getReflectionConditionFilter(ReflectionConditionFilter filter);

  std::vector<Kernel::V3D> getHKLs(double dMin, double dMax) const;
  std::vector<Kernel::V3D>
  getHKLs(double dMin, double dMax,
          HKLFilter_const_sptr reflectionConditionFilter) const;

  std::vector<Kernel::V3D> getUniqueHKLs(double dMin, double dMax) const;
  std::vector<Kernel::V3D>
  getUniqueHKLs(double dMin, double dMax,
                HKLFilter_const_sptr reflectionConditionFilter) const;

  std::vector<double> getDValues(const std::vector<Kernel::V3D> &hkls) const;
  std::vector<double> getFsSquared(const std::vector<Kernel::V3D> &hkls) const;

private:
  CrystalStructure m_crystalStructure;
  StructureFactorCalculator_sptr m_sfCalculator;
  HKLFilter_const_sptr m_defaultHKLFilter;
};

} // namespace Geometry
} // namespace Mantid

#endif /* MANTID_GEOMETRY_REFLECTIONGENERATOR_H_ */