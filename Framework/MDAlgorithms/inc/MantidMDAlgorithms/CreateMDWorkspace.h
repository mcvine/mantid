#ifndef MANTID_MDALGORITHMS_CREATEMDWORKSPACE_H_
#define MANTID_MDALGORITHMS_CREATEMDWORKSPACE_H_

#include "MantidAPI/Algorithm.h"
#include "MantidMDAlgorithms/BoxControllerSettingsAlgorithm.h"
#include "MantidDataObjects/MDEventFactory.h"
#include "MantidDataObjects/MDEventWorkspace.h"
#include "MantidGeometry/MDGeometry/MDFrame.h"

namespace Mantid {
namespace MDAlgorithms {

/** CreateMDWorkspace :
 *
 * Algorithm to create an empty MDEventWorkspace with a given number of
 *dimensions.
 *
 *
 * @author Janik Zikovsky
 * @date 2011-02-25 11:54:52.003137
 */
class DLLExport CreateMDWorkspace
    : public MDAlgorithms::BoxControllerSettingsAlgorithm {
public:
  CreateMDWorkspace();
  ~CreateMDWorkspace() override;

  /// Algorithm's name for identification
  const std::string name() const override { return "CreateMDWorkspace"; }
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Creates an empty MDEventWorkspace with a given number of "
           "dimensions.";
  }

  /// Algorithm's version for identification
  int version() const override { return 1; }
  /// Algorithm's category for identification
  const std::string category() const override {
    return "MDAlgorithms\\Creation";
  }
  std::map<std::string, std::string> validateInputs() override;

private:
  void init() override;
  void exec() override;

  template <typename MDE, size_t nd>
  void finish(typename DataObjects::MDEventWorkspace<MDE, nd>::sptr ws);
  Mantid::Geometry::MDFrame_uptr createMDFrame(std::string frame,
                                               std::string unit);
  bool checkIfFrameValid(const std::string &frame,
                         const std::vector<std::string> &targetFrames);
};

} // namespace Mantid
} // namespace DataObjects

#endif /* MANTID_MDALGORITHMS_CREATEMDWORKSPACE_H_ */
