#ifndef MANTID_ALGORITHMS_FILTEREVENTS_H_
#define MANTID_ALGORITHMS_FILTEREVENTS_H_

#include "MantidKernel/System.h"
#include "MantidAPI/Algorithm.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/SplittersWorkspace.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidAPI/ISplittersWorkspace.h"
#include "MantidKernel/TimeSplitter.h"
#include "MantidAPI/ITableWorkspace_fwd.h"
#include "MantidKernel/TimeSeriesProperty.h"

namespace Mantid {
namespace Algorithms {

class TimeAtSampleStrategy;

/** FilterEvents : Filter Events in EventWorkspace to multiple EventsWorkspace
  by Splitters

  @date 2012-04-04

  Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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
class DLLExport FilterEvents : public API::Algorithm {

  enum TOFCorrectionType {
    NoneCorrect,
    CustomizedCorrect,
    DirectCorrect,
    ElasticCorrect,
    IndirectCorrect
  };
  enum TOFCorrectionOp { MultiplyOp, ShiftOp };
  enum EVENTFILTERSKIP { EventFilterSkipNoDet, EventFilterSkipNoDetTOFCorr };

public:
  FilterEvents();
  ~FilterEvents() override;

  /// Algorithm's name for identification overriding a virtual method
  const std::string name() const override { return "FilterEvents"; }
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Filter events from an EventWorkspace to one or multiple "
           "EventWorkspaces according to a series of splitters.";
  }

  /// Algorithm's version for identification overriding a virtual method
  int version() const override { return 1; }

  /// Algorithm's category for identification overriding a virtual method
  const std::string category() const override {
    return "Events\\EventFiltering";
  }

private:
  // Implement abstract Algorithm methods
  void init() override;
  // Implement abstract Algorithm methods
  void exec() override;

  /// Process user input properties
  void processAlgorithmProperties();

  void processSplittersWorkspace();

  ///
  void processMatrixSplitterWorkspace();

  void createOutputWorkspaces();

  /// Set up detector calibration parameters
  void setupDetectorTOFCalibration();

  /// Set up detector calibration parameters for elastic scattering instrument
  TimeAtSampleStrategy *setupElasticTOFCorrection() const;

  /// Set up detector calibration parmaeters for direct inelastic scattering
  /// instrument
  TimeAtSampleStrategy *setupDirectTOFCorrection() const;

  /// Set up detector calibration parameters for indirect inelastic scattering
  /// instrument
  TimeAtSampleStrategy *setupIndirectTOFCorrection() const;

  /// Set up detector calibration parameters from customized values
  void setupCustomizedTOFCorrection();

  /// Filter events by splitters in format of Splitter
  void filterEventsBySplitters(double progressamount);

  /// Filter events by splitters in format of vector
  void filterEventsByVectorSplitters(double progressamount);

  /// Examine workspace
  void examineEventWS();

  DataObjects::EventWorkspace_sptr m_eventWS;
  DataObjects::SplittersWorkspace_sptr m_splittersWorkspace;
  API::MatrixWorkspace_const_sptr m_matrixSplitterWS;
  DataObjects::TableWorkspace_sptr m_detCorrectWorkspace;

  /// Flag to use matrix splitters or table splitters
  bool m_useTableSplitters;

  std::unordered_set<int> m_workGroupIndexes;
  Kernel::TimeSplitterType m_splitters;
  std::map<int, DataObjects::EventWorkspace_sptr> m_outputWS;
  std::vector<std::string> m_wsNames;

  std::vector<double> m_detTofOffsets;
  std::vector<double> m_detTofFactors;

  bool m_FilterByPulseTime;

  DataObjects::TableWorkspace_sptr m_informationWS;
  bool m_hasInfoWS;

  double m_progress;

  std::vector<std::string> getTimeSeriesLogNames();

  Kernel::TimeSplitterType generateSplitters(int wsindex);

  void splitLog(DataObjects::EventWorkspace_sptr eventws, std::string logname,
                Kernel::TimeSplitterType &splitters);

  /// Base of output workspace's name
  std::string m_outputWSNameBase;

  /// Flag to group workspace
  bool m_toGroupWS;

  /// Vector for splitting time
  std::vector<int64_t> m_vecSplitterTime;
  /// Vector for splitting grouip
  std::vector<int> m_vecSplitterGroup;

  /// Flag to split sample logs
  bool m_splitSampleLogs;

  /// Debug
  bool m_useDBSpectrum;
  int m_dbWSIndex;

  /// TOF detector/sample correction type
  TOFCorrectionType m_tofCorrType;

  /// Spectrum skip type
  EVENTFILTERSKIP m_specSkipType;
  /// Vector for skip information
  std::vector<bool> m_vecSkip;

  // Flag to have relative time in splitters workspace
  bool m_isSplittersRelativeTime;
  // Starting time for starting time of event filters
  Kernel::DateAndTime m_filterStartTime;
};

} // namespace Algorithms
} // namespace Mantid

#endif /* MANTID_ALGORITHMS_FILTEREVENTS_H_ */
