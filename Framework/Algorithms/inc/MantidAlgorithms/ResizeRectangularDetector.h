#ifndef MANTID_ALGORITHMS_RESIZERECTANGULARDETECTOR_H_
#define MANTID_ALGORITHMS_RESIZERECTANGULARDETECTOR_H_

#include "MantidKernel/System.h"
#include "MantidAPI/Algorithm.h"

namespace Mantid {
namespace Algorithms {

/** ResizeRectangularDetector : TODO: DESCRIPTION

  @date 2011-11-22

  Copyright &copy; 2011 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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
class DLLExport ResizeRectangularDetector : public API::Algorithm {
public:
  ResizeRectangularDetector();
  ~ResizeRectangularDetector() override;

  const std::string name() const override;
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Resize a RectangularDetector in X and/or Y.";
  }

  int version() const override;
  const std::string category() const override;

private:
  void init() override;
  void exec() override;
};

} // namespace Algorithms
} // namespace Mantid

#endif /* MANTID_ALGORITHMS_RESIZERECTANGULARDETECTOR_H_ */