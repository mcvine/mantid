#include "MantidWorkflowAlgorithms/LoadEventAndCompress.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/VisibleWhenProperty.h"

namespace Mantid {
namespace WorkflowAlgorithms {

using std::size_t;
using std::string;
using namespace Kernel;
using namespace API;
using namespace DataObjects;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(LoadEventAndCompress)

//----------------------------------------------------------------------------------------------
/** Constructor
 */
LoadEventAndCompress::LoadEventAndCompress() : m_filterBadPulses(EMPTY_DBL()) {}

//----------------------------------------------------------------------------------------------
/** Destructor
 */
LoadEventAndCompress::~LoadEventAndCompress() {}

//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const string LoadEventAndCompress::name() const { return "LoadEventAndCompress"; }

/// Algorithm's version for identification. @see Algorithm::version
int LoadEventAndCompress::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const string LoadEventAndCompress::category() const {
  return "Workflow\\DataHandling";
}

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const string LoadEventAndCompress::summary() const {
  return "Load an event workspace by chunks and compress";
}

//----------------------------------------------------------------------------------------------
/** Initialize the algorithm's properties.
 */
void LoadEventAndCompress::init() {
  // algorithms to copy properties from
    auto algLoadEventNexus = AlgorithmManager::Instance().createUnmanaged("LoadEventNexus");
    algLoadEventNexus->initialize();
    auto algDetermineChunking =
        AlgorithmManager::Instance().createUnmanaged("DetermineChunking");
    algDetermineChunking->initialize();

    // declare properties
    copyProperty(algLoadEventNexus, "Filename");
    copyProperty(algLoadEventNexus, "OutputWorkspace");
    copyProperty(algDetermineChunking, "MaxChunkSize");

    copyProperty(algLoadEventNexus, "FilterByTofMin");
    copyProperty(algLoadEventNexus, "FilterByTofMax");
    copyProperty(algLoadEventNexus, "FilterByTimeStart");
    copyProperty(algLoadEventNexus, "FilterByTimeStop");

    std::string grp1 = "Filter Events";
    setPropertyGroup("FilterByTofMin", grp1);
    setPropertyGroup("FilterByTofMax", grp1);
    setPropertyGroup("FilterByTimeStart", grp1);
    setPropertyGroup("FilterByTimeStop", grp1);

    copyProperty(algLoadEventNexus, "NXentryName");
    copyProperty(algLoadEventNexus, "LoadMonitors");
    copyProperty(algLoadEventNexus, "MonitorsAsEvents");
    copyProperty(algLoadEventNexus, "FilterMonByTofMin");
    copyProperty(algLoadEventNexus, "FilterMonByTofMax");
    copyProperty(algLoadEventNexus, "FilterMonByTimeStart");
    copyProperty(algLoadEventNexus, "FilterMonByTimeStop");

    setPropertySettings(
        "MonitorsAsEvents",
        new VisibleWhenProperty("LoadMonitors", IS_EQUAL_TO, "1"));
    IPropertySettings *asEventsIsOn =
        new VisibleWhenProperty("MonitorsAsEvents", IS_EQUAL_TO, "1");
    setPropertySettings("FilterMonByTofMin", asEventsIsOn);
    setPropertySettings("FilterMonByTofMax", asEventsIsOn->clone());
    setPropertySettings("FilterMonByTimeStart", asEventsIsOn->clone());
    setPropertySettings("FilterMonByTimeStop", asEventsIsOn->clone());

    std::string grp4 = "Monitors";
    setPropertyGroup("LoadMonitors", grp4);
    setPropertyGroup("MonitorsAsEvents", grp4);
    setPropertyGroup("FilterMonByTofMin", grp4);
    setPropertyGroup("FilterMonByTofMax", grp4);
    setPropertyGroup("FilterMonByTimeStart", grp4);
    setPropertyGroup("FilterMonByTimeStop", grp4);

    auto range = boost::make_shared<BoundedValidator<double>>();
    range->setBounds(0., 100.);
    declareProperty("FilterBadPulses", 95., range);
}

/// @see DataProcessorAlgorithm::determineChunk(const std::string &)
ITableWorkspace_sptr
LoadEventAndCompress::determineChunk(const std::string &filename) {
  double maxChunkSize = getProperty("MaxChunkSize");

  auto alg = createChildAlgorithm("DetermineChunking");
  alg->setProperty("Filename", filename);
  alg->setProperty("MaxChunkSize", maxChunkSize);
  alg->executeAsChildAlg();
  ITableWorkspace_sptr chunkingTable = alg->getProperty("OutputWorkspace");

  if (chunkingTable->rowCount() > 1)
    g_log.information() << "Will load data in " << chunkingTable->rowCount()
                        << " chunks\n";
  else
    g_log.information("Not chunking");

  return chunkingTable;
}

/// @see DataProcessorAlgorithm::loadChunk(const size_t)
MatrixWorkspace_sptr LoadEventAndCompress::loadChunk(const size_t rowIndex) {
  g_log.debug() << "loadChunk(" << rowIndex << ")\n";

  double rowCount = static_cast<double>(m_chunkingTable->rowCount());
  double progStart = static_cast<double>(rowIndex) / rowCount;
  double progStop = static_cast<double>(rowIndex + 1) / rowCount;

  auto alg = createChildAlgorithm("LoadEventNexus", progStart, progStop, true);
  alg->setProperty<string>("Filename", getProperty("Filename"));
  alg->setProperty<double>("FilterByTofMin", getProperty("FilterByTofMin"));
  alg->setProperty<double>("FilterByTofMax", getProperty("FilterByTofMax"));
  alg->setProperty<double>("FilterByTimeStart",
                           getProperty("FilterByTimeStart"));
  alg->setProperty<double>("FilterByTimeStop", getProperty("FilterByTimeStop"));

  alg->setProperty<string>("NXentryName", getProperty("NXentryName"));
  alg->setProperty<bool>("LoadMonitors", getProperty("LoadMonitors"));
  alg->setProperty<bool>("MonitorsAsEvents", getProperty("MonitorsAsEvents"));
  alg->setProperty<double>("FilterMonByTofMin",
                           getProperty("FilterMonByTofMin"));
  alg->setProperty<double>("FilterMonByTofMax",
                           getProperty("FilterMonByTofMax"));
  alg->setProperty<double>("FilterMonByTimeStart",
                           getProperty("FilterMonByTimeStart"));
  alg->setProperty<double>("FilterMonByTimeStop",
                           getProperty("FilterMonByTimeStop"));

  // set chunking information
  if (rowCount > 0.) {
    const std::vector<string> COL_NAMES = m_chunkingTable->getColumnNames();
    for (auto name = COL_NAMES.begin(); name != COL_NAMES.end(); ++name) {
      alg->setProperty(*name, m_chunkingTable->getRef<int>(*name, rowIndex));
    }
  }

  alg->executeAsChildAlg();
  Workspace_sptr wksp = alg->getProperty("OutputWorkspace");
  return boost::dynamic_pointer_cast<MatrixWorkspace>(wksp);
}

/**
 * Process a chunk in-place
 *
 * @param wksp
 */
void LoadEventAndCompress::processChunk(API::MatrixWorkspace_sptr wksp) {
  EventWorkspace_sptr eventWS =
      boost::dynamic_pointer_cast<EventWorkspace>(wksp);

  if (m_filterBadPulses > 0.) {
    auto filterBadPulses = createChildAlgorithm("FilterBadPulses");
    filterBadPulses->setProperty("InputWorkspace", eventWS);
    filterBadPulses->setProperty("OutputWorkspace", eventWS);
    filterBadPulses->setProperty("LowerCutoff", m_filterBadPulses);
    filterBadPulses->executeAsChildAlg();
  }

  auto compressEvents = createChildAlgorithm("CompressEvents");
  compressEvents->setProperty("InputWorkspace", eventWS);
  compressEvents->setProperty("OutputWorkspace", eventWS);
  compressEvents->executeAsChildAlg();
}

//----------------------------------------------------------------------------------------------
/** Execute the algorithm.
 */
void LoadEventAndCompress::exec() {
  m_filename = getPropertyValue("Filename");
  m_filterBadPulses = getProperty("FilterBadPulses");

  m_chunkingTable = determineChunk(m_filename);

  // first run is free
  MatrixWorkspace_sptr resultWS = loadChunk(0);
  processChunk(resultWS);

  // load the other chunks
  const size_t numRows = m_chunkingTable->rowCount();
  for (size_t i = 1; i < numRows; ++i) {
    MatrixWorkspace_sptr temp = loadChunk(i);
    processChunk(temp);
    resultWS = plus2(resultWS,temp);
  }
  Workspace_sptr total = assemble(resultWS);

  // Don't bother compressing combined workspace. DetermineChunking is designed
  // to prefer loading full banks so no further savings should be available.

  setProperty("OutputWorkspace", total);
}

} // namespace WorkflowAlgorithms
} // namespace Mantid
