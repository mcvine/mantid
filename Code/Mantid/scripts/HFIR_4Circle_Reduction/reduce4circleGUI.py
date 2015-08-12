#pylint: disable=invalid-name
################################################################################
#
# MainWindow application for reducing HFIR 4-circle 
#
################################################################################
import sys
import os
import csv

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import reduce4circleControl as r4c


try:
    import mantid
except ImportError:
    sys.path.append('/home/wzz/Mantid/Code/debug/bin/')
    import mantid
finally:
    import mantid.simpleapi as api
    import mantid.kernel
    from mantid.simpleapi import AnalysisDataService
    from mantid.kernel import ConfigService

import fourcircle_utility as fcutil


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

from Ui_MainWindow import Ui_MainWindow #import line for the UI python class


class MainWindow(QtGui.QMainWindow):
    """ Class of Main Window (top)
    """
    def __init__(self, parent=None):
        """ Initialization and set up
        """
        # Base class
        QtGui.QMainWindow.__init__(self,parent)

        # UI Window (from Qt Designer)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Mantid configuration
        self._instrument = str(self.ui.comboBox_instrument.currentText())
        # config = ConfigService.Instance()
        # self._instrument = config["default.instrument"]

        # Event handling definitions
        # Tab 'Data Access'
        self.connect(self.ui.pushButton_browseLocalDataDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_spice_data)
        self.connect(self.ui.pushButton_testURLs, QtCore.SIGNAL('clicked()'),
                     self.do_test_url)
        self.connect(self.ui.pushButton_ListScans, QtCore.SIGNAL('clicked()'),
                     self.do_list_scans)
        self.connect(self.ui.pushButton_downloadExpData, QtCore.SIGNAL('clicked()'),
                     self.do_download_spice_data)
        self.connect(self.ui.comboBox_mode, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_data_access_mode)

        # Tab 'View Raw Data'
        self.connect(self.ui.pushButton_plotRawPt, QtCore.SIGNAL('clicked()'),
                     self.do_plot_pt_raw)
        self.connect(self.ui.pushButton_prevPtNumber, QtCore.SIGNAL('clicked()'),
                     self.do_plot_prev_pt_raw)
        self.connect(self.ui.pushButton_nextPtNumber, QtCore.SIGNAL('clicked()'),
                     self.do_plot_next_pt_raw)
        self.connect(self.ui.pushButton_showPtList, QtCore.SIGNAL('clicked()'),
                     self.show_scan_pt_list)

        # Tab 'Advanced'
        self.connect(self.ui.pushButton_browseLocalCache, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_cache_dir)
        self.connect(self.ui.pushButton_browseWorkDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_working_dir)
        self.connect(self.ui.comboBox_instrument, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_instrument_name)

        # Menu
        self.connect(self.ui.actionSave_Session, QtCore.SIGNAL('triggered()'),
                     self.save_current_session)
        self.connect(self.ui.actionLoad_Session, QtCore.SIGNAL('triggered()'),
                     self.load_session)

        # Tab 'calculate ub matrix'
        self.connect(self.ui.pushButton_findPeak, QtCore.SIGNAL('clicked()'),
                     self.do_find_peak)

        self.connect(self.ui.pushButton_calUB, QtCore.SIGNAL('clicked()'),
                self.doCalUBMatrix)

        self.connect(self.ui.pushButton_acceptUB, QtCore.SIGNAL('clicked()'),
                self.doAcceptCalUB)

        self.connect(self.ui.pushButton_resetCalUB, QtCore.SIGNAL('clicked()'),
                self.doResetCalUB) 

        # Event handling for tab 'refine ub matrix'
        self.connect(self.ui.pushButton_addToRefine, QtCore.SIGNAL('clicked()'),
                self.doAddScanPtToRefineUB)

        # Validator


        # Declaration of class variable
        self._runID = None
        self._expID = None
        self._currPt = None
        self._xmlwsbasename = None
        
        # Some configuration
        self._homeSrcDir = os.getcwd()
        self._homeDir = os.getcwd()

        # Control
        self._myControl = r4c.CWSCDReductionControl(self._instrument)
        self._allowDownload = True
        self._dataAccessMode = 'Download'

        # Initial setup
        self.ui.tabWidget.setCurrentIndex(0)

        # Tab 'Access'
        self.ui.lineEdit_url.setText('http://neutron.ornl.gov/user_data/hb3a/')
        self.ui.comboBox_mode.setCurrentIndex(0)
        self.ui.lineEdit_localSpiceDir.setEnabled(False)
        self.ui.pushButton_browseLocalDataDir.setEnabled(False)

        return

    def change_data_access_mode(self):
        """ Change data access mode between downloading from server and local
        Event handling methods
        :return:
        """
        new_mode = str(self.ui.comboBox_mode.currentText())
        self._dataAccessMode = new_mode

        if new_mode.startswith('Local') is True:
            self.ui.lineEdit_localSpiceDir.setEnabled(True)
            self.ui.pushButton_browseLocalDataDir.setEnabled(True)
            self.ui.lineEdit_url.setEnabled(False)
            self.ui.lineEdit_localSrcDir.setEnabled(False)
            self.ui.pushButton_browseLocalCache.setEnabled(False)
            self._allowDownload = False
        else:
            self.ui.lineEdit_localSpiceDir.setEnabled(False)
            self.ui.pushButton_browseLocalDataDir.setEnabled(False)
            self.ui.lineEdit_url.setEnabled(True)
            self.ui.lineEdit_localSrcDir.setEnabled(True)
            self.ui.pushButton_browseLocalCache.setEnabled(True)
            self._allowDownload = True

        return

    def change_instrument_name(self):
        """ Handing the event as the instrument name is changed
        :return:
        """
        new_instrument = str(self.ui.comboBox_instrument.currentText())
        self.pop_one_button_dialog('Change of instrument during data processing is dangerous.')
        status, error_message = self._myControl.set_instrument_name(new_instrument)
        if status is False:
            self.pop_one_button_dialog(error_message)

        return

    def doAcceptCalUB(self):
        """ Accept the calculated UB matrix
        """

        return

    def doAddScanPtToRefineUB(self):
        """ Add scan/pt numbers to the list of data points for refining ub matrix

        And the added scan number and pt numbers will be reflected in the (left sidebar)

        """
        raise NotImplementedError("ASAP")

        return

    def do_browse_local_cache_dir(self):
        """ Browse local cache directory
        :return:
        """
        local_cache_dir = str(QtGui.QFileDialog.getExistingDirectory(self,
                                                                     'Get Local Cache Directory',
                                                                     self._homeSrcDir))

        # Set local directory to control
        status, error_message = self._myControl.set_local_data_dir(local_cache_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        # Synchronize to local data/spice directory and local cache directory
        if str(self.ui.lineEdit_localSpiceDir.text()) != '':
            prev_dir = str(self.ui.lineEdit_localSrcDir.text())
            self.pop_one_button_dialog('Local data directory was set up as %s' %
                                       prev_dir)
        self.ui.lineEdit_localSrcDir.setText(local_cache_dir)
        self.ui.lineEdit_localSpiceDir.setText(local_cache_dir)

        return

    def do_browse_local_spice_data(self):
        """ Browse local source SPICE data directory
        """
        src_spice_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Directory',
                                                                   self._homeSrcDir))
        # Set local data directory to controller
        status, error_message = self._myControl.set_local_data_dir(src_spice_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        self._homeSrcDir = src_spice_dir
        self.ui.lineEdit_localSpiceDir.setText(src_spice_dir)

        return

    def do_browse_working_dir(self):
        """
        Browse and set up working directory
        :return:
        """
        work_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Working Directory'), self._homeDir)
        status, error_message = self._myControl.set_working_directory(work_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)

        return
    
    def doCalUBMatrix(self):
        """ Calculate UB matrix by 2 or 3 reflections
        """

        return

    def do_download_spice_data(self):
        """ Download SPICE data
        :return:
        """
        # Check scans to download
        scan_list_str = str(self.ui.lineEdit_downloadScans.text())
        if len(scan_list_str) > 0:
            # user specifies scans to download
            valid, scan_list = fcutil.parse_int_array(scan_list_str)
            if valid is False:
                error_message = scan_list
                self.pop_one_button_dialog(scan_list)
        else:
            # Get all scans
            # FIXME - Implement all scan case
            scan_list = fcutil.get_scans_list(server_url, exp_no, return_list=True)
        self.pop_one_button_dialog('Going to download scans %s.' % str(scan_list))

        # Check location
        destination_dir = str(self.ui.lineEdit_localSrcDir.text())
        status, error_message = self._myControl.set_local_data_dir(destination_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
        else:
            self.pop_one_button_dialog('Spice files will be downloaded to %s.' % destination_dir)

        # Set up myControl for downloading data
        exp_no = int(self.ui.lineEdit_exp.text())
        self._myControl.set_exp_number(exp_no)

        server_url = str(self.ui.lineEdit_url.text())
        status, error_message = self._myControl.set_server_url(server_url)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        # Download
        self._myControl.download_data_set(scan_list)

        return
    
    def do_find_peak(self):
        """ Find peak in a given scan/pt and record it
        """
        # Get experiment, scan and pt
        status, ret_obj = self._parse_integers_editor([self.ui.lineEdit_exp, self.ui.lineEdit_scanNumber,
                                                       self.ui.lineEdit_ptNumber])
        if status is True:
            exp_no, scan_no, pt_no = ret_obj
        else:
            self.pop_one_button_dialog(ret_obj)

        # Find peak
        status, error_message = self._myProject.findPeak(exp_no, scan_no, pt_no)
        if status is False:
            self.pop_one_button_dialog(error_message)

        # Set up correct values

        return

    def do_list_scans(self):
        """ List all scans available
        :return:
        """
        # Experiment number
        exp_no = int(self.ui.lineEdit_exp.text())

        access_mode = str(self.ui.comboBox_mode.currentText())
        if access_mode == 'Local':
            spice_dir = str(self.ui.lineEdit_localSpiceDir.text())
            message = fcutil.get_scans_list_local_disk(spice_dir, exp_no)
        else:
            url = str(self.ui.lineEdit_url.text())
            message = fcutil.get_scans_list(url, exp_no)

        self.pop_one_button_dialog(message)

        return

    def do_plot_pt_raw(self):
        """ Plot the Pt. 
        """
        # Get measurement pt and the file number
        status, ret_obj = self._parse_integers_editor([self.ui.lineEdit_exp,
                                                       self.ui.lineEdit_run,
                                                       self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Call to plot 2D
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return
        
    def do_plot_prev_pt_raw(self):
        """ Plot the Pt. 
        """
        # Get measurement pt and the file number
        status, ret_obj = self._parse_integers_editor([self.ui.lineEdit_exp,
                                                       self.ui.lineEdit_run,
                                                       self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Previous one
        pt_no -= 1
        if pt_no <= 0:
            self.pop_one_button_dialog('Pt. = 1 is the first one.')
            return
        else:
            self.ui.lineEdit_rawDataPtNo.setText('%d' % pt_no)

        # Plot
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def do_plot_next_pt_raw(self):
        """ Plot the Pt. 
        """
        # Get measurement pt and the file number
        status, ret_obj = self._parse_integers_editor([self.ui.lineEdit_exp,
                                                       self.ui.lineEdit_run,
                                                       self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Previous one
        pt_no += 1
        # get last Pt. number
        status, last_pt_no = self._myControl.get_pt_numbers(exp_no, scan_no)[-1]
        if status is False:
            error_message = last_pt_no
            self.pop_one_button_dialog('Unable to access Spice table for scan %d due to %s.' % (scan_no, error_message))
        if pt_no > last_pt_no:
            self.pop_one_button_dialog('Pt. = %d is the last one of scan %d.' % (pt_no, scan_no))
            return
        else:
            self.ui.lineEdit_rawDataPtNo.setText('%d' % pt_no)

        # Plot
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def doResetCalUB(self):
        """ Reset/reject the UB matrix calculation
        """

        return

    def do_test_url(self):
        """ Test whether the root URL provided specified is good
        """
        url = str(self.ui.lineEdit_url.text())

        url_is_good = fcutil.check_url(url)
        if url_is_good is True:
            self.pop_one_button_dialog("URL %s is valid." % url)
        else:
            self.pop_one_button_dialog("Unable to access %s.  Check internet access." % url)

        return

    def pop_one_button_dialog(self, message):
        """ Pop up a one-button dialog
        :param message:
        :return:
        """
        QtGui.QMessageBox.information(self, '4-circle Data Reduction', message)

        return

    def save_current_session(self, filename=None):
        """ Save current session/value setup to
        :return:
        """
        # Set up dictionary
        save_dict = {}
        save_dict['lineEdit_exp'] = str(self.ui.lineEdit_exp.text())
        save_dict['lineEdit_localSpiceDir'] = str(self.ui.lineEdit_localSpiceDir.text())

        save_dict['comboBox_mode'] = self.ui.comboBox_mode.currentIndex()

        # Save to csv file
        if filename is None:
            filename = 'session_backup.csv'
        ofile = open(filename, 'w')
        writer = csv.writer(ofile)
        for key, value in save_dict.items():
            writer.writerow([key, value])
        ofile.close()

        return

    def load_session(self, filename=None):
        """
        To read it back:
        :param filename:
        :return:
        """
        # TODO - Doc!
        if filename is None:
            filename = 'session_backup.csv'

        in_file = open(filename, 'r')
        reader = csv.reader(in_file)
        my_dict = dict(x for x in reader)

        # ...
        for key, value in my_dict.items():
            if key.startswith('lineEdit') is True:
                self.ui.__getattribute__(key).setText(value)
            elif key.startswith('comboBox') is True:
                self.ui.__getattribute__(key).setCurrentIndex(int(value))
            else:
                self.pop_one_button_dialog('Error! Widget name %s is not supported' % key)
        # END-FOR

        # ...
        self._myControl.set_local_data_dir(str(self.ui.lineEdit_localSpiceDir.text()))

        return

    def show_scan_pt_list(self):
        """ Show the range of Pt. in a scan
        :return:
        """
        # Get parameters
        status, inp_list = self._parse_integers_editor([self.ui.lineEdit_exp, self.ui.lineEdit_run])
        if status is False:
            self.pop_one_button_dialog(inp_list)
            return
        else:
            exp_no = inp_list[0]
            scan_no = inp_list[1]

        status, ret_obj = self._myControl.get_pt_numbers(exp_no, scan_no)

        # Form message
        if status is False:
            # Failed to get Pt. list
            error_message = ret_obj
            self.pop_one_button_dialog(error_message)
        else:
            # Form message
            pt_list = sorted(ret_obj)
            num_pts = len(pt_list)
            info = 'Exp %d Scan %d has %d Pt. ranging from %d to %d.\n' % (exp_no, scan_no, num_pts,
                                                                           pt_list[0], pt_list[-1])
            num_miss_pt = pt_list[-1] - pt_list[0] + 1 - num_pts
            if num_miss_pt > 0:
                info += 'There are %d Pt. skipped.\n' % num_miss_pt

            self.pop_one_button_dialog(info)

        return

    def _plot_raw_xml_2d(self, exp_no, scan_no, pt_no):
        """ Plot raw workspace from XML file for a measurement/pt.
        """
        # Load Data including SPICE scan file (necessary???) and Pt's xml file
        status, error_message = self._myControl.load_spice_scan_file(exp_no, scan_no)
        if status is False and self._allowDownload is False:
            self.pop_one_button_dialog(error_message)
            return
        elif status is False and self._allowDownload is True:
            status, error_message = self._myControl.download_spice_file(exp_no, scan_no)
            if status is False:
                self.pop_one_button_dialog(error_message)
                return

        status, error_message = self._myControl.load_spice_xml_file(exp_no, scan_no, pt_no)
        if status is False:
            if self._allowDownload is True:
                status, error_message = self._myControl.download_spice_xml_file(exp_no, scan_no)
                if status is False:
                    self.pop_one_button_dialog(error_message)
                    return
            else:
                self.pop_one_button_dialog(error_message)

        # Convert a list of vector to 2D numpy array for imshow()
        # Get data and plot
        raw_det_data = self._myControl.get_raw_detector_counts(exp_no, scan_no, pt_no)
        self.ui.graphicsView.clear_canvas()
        self.ui.graphicsView.add_plot_2d(raw_det_data, x_min=0, x_max=256, y_min=0, y_max=256,
                                         hold_prev_image=False)

        return

    def _parse_integers_editor(self, line_edit_list):
        """
        :param line_edit_list:
        :return: (True, list of integers); (False, error message)
        """
        error_message = ''
        integer_list = []

        for line_edit in line_edit_list:
            try:
                str_value = str(line_edit.text()).strip()
                int_value = int(str_value)
            except ValueError as e:
                error_message += 'Unable to parse to integer. %s\n' % (str(e))
            else:
                if str_value != '%d' % int_value:
                    error_message += 'Value %s is not a proper integer.\n' % str_value
                else:
                    integer_list.append(int_value)
                    print 'Value %s to %d' % (str_value, int_value)
            # END-TRY
        # END-FOR

        if len(error_message) > 0:
            self.pop_one_button_dialog(error_message)
            return False, error_message

        return True, integer_list
