//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAlgorithms/ConvertToDistribution.h"
#include "MantidAPI/WorkspaceValidators.h"

namespace Mantid
{
namespace Algorithms
{

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(ConvertToDistribution)

using namespace API;

void ConvertToDistribution::init()
{
  CompositeValidator<> *wsValidator = new CompositeValidator<>;
  wsValidator->add(new HistogramValidator<>);
  wsValidator->add(new RawCountValidator<>);
  declareProperty(new WorkspaceProperty<>("Workspace","",Kernel::Direction::InOut,wsValidator));
}

void ConvertToDistribution::exec()
{
  Workspace_sptr workspace = getProperty("Workspace");

  WorkspaceHelpers::makeDistribution(workspace);
}

} // namespace Algorithms
} // namespace Mantid
