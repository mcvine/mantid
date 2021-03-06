#pylint: disable=no-init
import mantid.simpleapi as api
from mantid.api import *
from mantid.kernel import *
from mantid import config

MICROEV_TO_MILLIEV = 1000.0
DEFAULT_BINS = [0., 0., 0.]
DEFAULT_RANGE = [6.24, 6.30]
DEFAULT_MASK_GROUP_DIR = "/SNS/BSS/shared/autoreduce"
DEFAULT_MASK_FILE = "BASIS_Mask.xml"
DEFAULT_VANADIUM_ENERGY_RANGE = [-0.0034, 0.0034]  # meV

#pylint: disable=too-many-instance-attributes
class BASISReduction(PythonAlgorithm):

    _short_inst = None
    _long_inst = None
    _extension = None
    _doIndiv = None
    _etBins = None
    _qBins = None
    _noMonNorm = None
    _maskFile = None
    _groupDetOpt = None
    _overrideMask = None
    _dMask = None

    # class variables related to division by Vanadium (normalization)
    _normRange = None
    _norm_run_list = None
    _normWs = None
    _normMonWs = None

    _run_list = None  # a list of runs, or a list of sets of runs
    _samWs = None
    _samMonWs = None
    _samWsRun = None
    _samSqwWs = None

    def __init__(self):
        PythonAlgorithm.__init__(self)
        self._doNorm = None  # stores the selected item from normalization_types
        self._normalizeToFirst = False

    def category(self):
        return "Inelastic\\Reduction"

    def name(self):
        return "BASISReduction"

    def version(self):
        return 1

    def summary(self):
        return "Multiple-file BASIS reduction."

    def PyInit(self):
        self._short_inst = "BSS"
        self._long_inst = "BASIS"
        self._extension = "_event.nxs"

        self.declareProperty("RunNumbers", "", "Sample run numbers")
        self.declareProperty("DoIndividual", False, "Do each run individually")
        self.declareProperty("NoMonitorNorm", False,
                             "Stop monitor normalization")
        arrVal = FloatArrayLengthValidator(2)

        self.declareProperty(FloatArrayProperty("EnergyBins", DEFAULT_BINS,
                                                direction=Direction.Input),
                             "Energy transfer binning scheme (in ueV)")
        self.declareProperty(FloatArrayProperty("MomentumTransferBins",
                                                DEFAULT_BINS,
                                                direction=Direction.Input),
                             "Momentum transfer binning scheme")
        self.declareProperty(FileProperty(name="MaskFile", defaultValue="",
                                          action=FileAction.OptionalLoad, extensions=['.xml']),
                             "Directory location for standard masking and grouping files.")
        grouping_type = ["None", "Low-Resolution", "By-Tube"]
        self.declareProperty("GroupDetectors", "None",
                             StringListValidator(grouping_type),
                             "Switch for grouping detectors")

        self.declareProperty("NormalizeToFirst", False, "Normalize spectra to intensity of spectrum with lowest Q?")

        # Properties setting the division by vanadium
        titleDivideByVanadium = "Normalization by Vanadium"
        self.declareProperty("DivideByVanadium", False, direction=Direction.Input,
                             doc="Do we normalize by the vanadium intensity?")
        self.setPropertyGroup("DivideByVanadium", titleDivideByVanadium)
        ifDivideByVanadium = EnabledWhenProperty("DivideByVanadium",
                                                 PropertyCriterion.IsNotDefault)

        normalization_types = ["by Q slice", "by detector ID"]
        self.declareProperty("NormalizationType", "by Q slice",
                             StringListValidator(normalization_types),
                             "Select a Vanadium normalization")
        self.setPropertySettings("NormalizationType", ifDivideByVanadium)
        self.setPropertyGroup("NormalizationType", titleDivideByVanadium)

        self.declareProperty("NormRunNumbers", "", "Normalization run numbers")
        self.setPropertySettings("NormRunNumbers", ifDivideByVanadium)
        self.setPropertyGroup("NormRunNumbers", titleDivideByVanadium)

        self.declareProperty(FloatArrayProperty("NormWavelengthRange", DEFAULT_RANGE,
                                                arrVal, direction=Direction.Input),
                             "Wavelength range for normalization")
        self.setPropertySettings("NormWavelengthRange", ifDivideByVanadium)
        self.setPropertyGroup("NormWavelengthRange", titleDivideByVanadium)

    def PyExec(self):
        config['default.facility'] = "SNS"
        config['default.instrument'] = self._long_inst
        self._doIndiv = self.getProperty("DoIndividual").value
        self._etBins = self.getProperty("EnergyBins").value / MICROEV_TO_MILLIEV
        self._qBins = self.getProperty("MomentumTransferBins").value
        self._noMonNorm = self.getProperty("NoMonitorNorm").value
        self._maskFile = self.getProperty("MaskFile").value
        self._groupDetOpt = self.getProperty("GroupDetectors").value
        self._normalizeToFirst = self.getProperty("NormalizeToFirst").value
        self._normalizeToVanadium = self.getProperty("GroupDetectors").value
        self._doNorm = self.getProperty("DivideByVanadium").value

        datasearch = config["datasearch.searcharchive"]
        if datasearch != "On":
            config["datasearch.searcharchive"] = "On"

        # Handle masking file override if necessary
        self._overrideMask = bool(self._maskFile)
        if not self._overrideMask:
            config.appendDataSearchDir(DEFAULT_MASK_GROUP_DIR)
            self._maskFile = DEFAULT_MASK_FILE

        api.LoadMask(Instrument='BASIS', OutputWorkspace='BASIS_MASK',
                     InputFile=self._maskFile)

        # Work around length issue
        _dMask = api.ExtractMask('BASIS_MASK')
        self._dMask = _dMask[1]
        api.DeleteWorkspace(_dMask[0])

        ############################
        ##  Process the Vanadium  ##
        ############################

        norm_runs = self.getProperty("NormRunNumbers").value
        if self._doNorm and bool(norm_runs):
            if ";" in norm_runs:
                raise SyntaxError("Normalization does not support run groups")
            self._doNorm = self.getProperty("NormalizationType").value
            self.log().information("Divide by Vanadium with normalization" + self._doNorm)

            # The following steps are common to all types of Vanadium normalization

            # norm_runs encompasses a single set, thus _getRuns returns
            # a list of only one item
            norm_set = self._getRuns(norm_runs, doIndiv=False)[0]
            self._normWs = self._sum_and_calibrate(norm_set, extra_extension="_norm")

            # This rebin integrates counts onto a histogram of a single bin
            if self._doNorm == "by detectorID":
                normRange = self.getProperty("NormWavelengthRange").value
                self._normRange = [normRange[0], normRange[1]-normRange[0], normRange[1]]
                api.Rebin(InputWorkspace=self._normWs, OutputWorkspace=self._normWs,
                          Params=self._normRange)

            # FindDetectorsOutsideLimits to be substituted by MedianDetectorTest
            api.FindDetectorsOutsideLimits(InputWorkspace=self._normWs,
                                           OutputWorkspace="BASIS_NORM_MASK")

            # additional reduction steps when normalizing by Q slice
            if self._doNorm == "by Q slice":
                self._normWs = self._group_and_SofQW(self._normWs, self._etBins, isSample=False)

        ##########################
        ##  Process the sample  ##
        ##########################
        self._run_list = self._getRuns(self.getProperty("RunNumbers").value,
                                       doIndiv=self._doIndiv)
        for run_set in self._run_list:
            self._samWs = self._sum_and_calibrate(run_set)
            self._samWsRun = str(run_set[0])
            # Mask detectors with insufficient Vanadium signal
            if self._doNorm:
                api.MaskDetectors(Workspace=self._samWs,
                                  MaskedWorkspace='BASIS_NORM_MASK')
            # Divide by Vanadium
            if self._doNorm == "by detector ID":
                api.Divide(LHSWorkspace=self._samWs, RHSWorkspace=self._normWs,
                           OutputWorkspace=self._samWs)
            # additional reduction steps
            self._samSqwWs = self._group_and_SofQW(self._samWs, self._etBins, isSample=True)
            # Divide by Vanadium
            if self._doNorm == "by Q slice":
                api.Integration(InputWorkspace=self._normWs, OutputWorkspace=self._normWs,
                                RangeLower=DEFAULT_VANADIUM_ENERGY_RANGE[0],
                                RangeUpper=DEFAULT_VANADIUM_ENERGY_RANGE[1])
                api.Divide(LHSWorkspace=self._samSqwWs, RHSWorkspace=self._normWs,
                           OutputWorkspace=self._samSqwWs)
            # Clear mask from reduced file. Needed for binary operations
            # involving this S(Q,w)
            api.ClearMaskFlag(Workspace=self._samSqwWs)
            # Scale so that elastic line has Y-values ~ 1
            if self._normalizeToFirst:
                self._ScaleY(self._samSqwWs)
            # Output Dave and Nexus files
            extension = "_divided.dat" if self._doNorm else ".dat"
            dave_grp_filename = self._makeRunName(self._samWsRun,
                                                  False) + extension
            api.SaveDaveGrp(Filename=dave_grp_filename,
                            InputWorkspace=self._samSqwWs,
                            ToMicroEV=True)
            extension = "_divided_sqw.nxs" if self._doNorm else "_sqw.nxs"
            processed_filename = self._makeRunName(self._samWsRun,
                                                   False) + extension
            api.SaveNexus(Filename=processed_filename,
                          InputWorkspace=self._samSqwWs)

    def _getRuns(self, rlist, doIndiv=True):
        """
        Create sets of run numbers for analysis. A semicolon indicates a
        separate group of runs to be processed together.
        @param rlist: string containing all the run numbers to be reduced.
        @return if _doIndiv is False, return a list of IntArrayProperty objects.
         Each item is a pseudolist containing a set of runs to be reduced together.
         if _doIndiv is True, return a list of strings, each string is a run number.
        """
        run_list = []
        # ";" separates the runs into substrings. Each substring represents a set of runs
        rlvals = rlist.split(';')
        for rlval in rlvals:
            iap = IntArrayProperty("", rlval)  # split the substring
            if doIndiv:
                run_list.extend([[x] for x in iap.value])
            else:
                run_list.append(iap.value)
        return run_list

    def _makeRunName(self, run, useShort=True):
        """
        Make name like BSS_24234
        """
        if useShort:
            return self._short_inst + "_" + str(run)
        else:
            return self._long_inst + "_" + str(run)

    def _makeRunFile(self, run):
        """
        Make name like BSS24234
        """
        return self._short_inst + str(run)

    def _sumRuns(self, run_set, sam_ws, mon_ws, extra_ext=None):
        """
        Aggregate the set of runs
        @param run_set: list of run numbers
        @param sam_ws:  name of aggregate workspace for the sample
        @param mon_ws:  name of aggregate workspace for the monitors
        @param extra_ext: string to be added to the temporary workspaces
        """
        for run in run_set:
            ws_name = self._makeRunName(run)
            if extra_ext is not None:
                ws_name += extra_ext
            mon_ws_name = ws_name  + "_monitors"
            run_file = self._makeRunFile(run)

            api.Load(Filename=run_file, OutputWorkspace=ws_name)
            if not self._noMonNorm:
                api.LoadNexusMonitors(Filename=run_file,
                                      OutputWorkspace=mon_ws_name)
            if sam_ws != ws_name:
                api.Plus(LHSWorkspace=sam_ws, RHSWorkspace=ws_name,
                         OutputWorkspace=sam_ws)
                api.DeleteWorkspace(ws_name)
            if mon_ws != mon_ws_name and not self._noMonNorm:
                api.Plus(LHSWorkspace=mon_ws,
                         RHSWorkspace=mon_ws_name,
                         OutputWorkspace=mon_ws)
                api.DeleteWorkspace(mon_ws_name)

    def _calibData(self, sam_ws, mon_ws):
        api.MaskDetectors(Workspace=sam_ws,
                          DetectorList=self._dMask)
                          #MaskedWorkspace='BASIS_MASK')
        api.ModeratorTzeroLinear(InputWorkspace=sam_ws,\
                           OutputWorkspace=sam_ws)
        api.LoadParameterFile(Workspace=sam_ws,
                              Filename=config.getInstrumentDirectory() + 'BASIS_silicon_111_Parameters.xml')
        api.ConvertUnits(InputWorkspace=sam_ws,
                         OutputWorkspace=sam_ws,
                         Target='Wavelength', EMode='Indirect')

        if not self._noMonNorm:
            api.ModeratorTzeroLinear(InputWorkspace=mon_ws,\
                               OutputWorkspace=mon_ws)
            api.Rebin(InputWorkspace=mon_ws,
                      OutputWorkspace=mon_ws, Params='10')
            api.ConvertUnits(InputWorkspace=mon_ws,
                             OutputWorkspace=mon_ws,
                             Target='Wavelength')
            api.OneMinusExponentialCor(InputWorkspace=mon_ws,
                                       OutputWorkspace=mon_ws,
                                       C='0.20749999999999999',
                                       C1='0.001276')
            api.Scale(InputWorkspace=mon_ws,
                      OutputWorkspace=mon_ws,
                      Factor='9.9999999999999995e-07')
            api.RebinToWorkspace(WorkspaceToRebin=sam_ws,
                                 WorkspaceToMatch=mon_ws,
                                 OutputWorkspace=sam_ws)
            api.Divide(LHSWorkspace=sam_ws,
                       RHSWorkspace=mon_ws,
                       OutputWorkspace=sam_ws)

    def _sum_and_calibrate(self, run_set, extra_extension=""):
        """
        Aggregate the set of runs and calibrate
        @param run_set: list of run numbers
        @param extra_extension: string to be added to the workspace names
        @return: workspace name of the aggregated and calibrated data
        """
        wsName = self._makeRunName(run_set[0])
        wsName += extra_extension
        wsMonName = wsName + "_monitors"
        self._sumRuns(run_set, wsName, wsMonName, extra_extension)
        self._calibData(wsName, wsMonName)
        return wsName

    def _group_and_SofQW(self, wsName, etRebins, isSample=True):
        """ Transforms from wavelength and detector ID to S(Q,E)
        @param wsName: workspace as a function of wavelength and detector id
        @return: S(Q,E)
        """
        api.ConvertUnits(InputWorkspace=wsName,
                         OutputWorkspace=wsName,
                         Target='DeltaE', EMode='Indirect')
        api.CorrectKiKf(InputWorkspace=wsName,
                        OutputWorkspace=wsName,
                        EMode='Indirect')
        api.Rebin(InputWorkspace=wsName,
                  OutputWorkspace=wsName,
                  Params=etRebins)
        if self._groupDetOpt != "None":
            if self._groupDetOpt == "Low-Resolution":
                grp_file = "BASIS_Grouping_LR.xml"
            else:
                grp_file = "BASIS_Grouping.xml"
            # If mask override used, we need to add default grouping file
            # location to search paths
            if self._overrideMask:
                config.appendDataSearchDir(DEFAULT_MASK_GROUP_DIR)
                api.GroupDetectors(InputWorkspace=wsName,
                                   OutputWorkspace=wsName,
                                   MapFile=grp_file, Behaviour="Sum")
        wsSqwName = wsName+'_divided_sqw' if isSample and self._doNorm else wsName+'_sqw'
        api.SofQW3(InputWorkspace=wsName,
                   OutputWorkspace=wsSqwName,
                   QAxisBinning=self._qBins, EMode='Indirect',
                   EFixed='2.0826')
        return wsSqwName

    def _ScaleY(self, wsName):
        """
        Scale all spectra by a number so that the maximum of the first spectra
        is rescaled to 1
        @param wsName: name of the workspace to rescale
        """
        workspace = mtd[wsName]
        maximumYvalue = workspace.dataY(0).max()
        api.Scale(InputWorkspace=wsName,
                  OutputWorkspace=wsName,
                  Factor=1./maximumYvalue, Operation="Multiply",)

# Register algorithm with Mantid.
AlgorithmFactory.subscribe(BASISReduction)
