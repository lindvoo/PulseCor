#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 15:24:09 2019

@author: lindvoo
"""

# General
import os
import sys
import numpy as np

# Basis for the GUI
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QFileDialog, QInputDialog, QLineEdit, QSlider, QLabel

# Used to plot the data
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.qt_compat import QtCore
from matplotlib.figure import Figure

# Peak detection
from scipy.ndimage.filters import maximum_filter

# Convert EEG files
import mne

# Convert ACQ files
import bioread

# Filter
from scipy.signal import butter, filtfilt

os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")

# ------------------------------------------------------------------------------

# Settings
class DATASETTINGS():
    def __init__(self):
        self.Hz = 50  # recording frequency [.eeg files and .acq files are downsampled to 50 hz]
        self.baseline = 1  # seconds > can be changed in GUI
        self.stimduur = 6  # seconds > can be changed in GUI


# ------------------------------------------------------------------------------

# Main window
class Window(QMainWindow):

    def __init__(self):

        super().__init__()

        title = "PulseCor"
        top = 0
        left = 0
        self.width = 1300
        self.height = 600

        self.setWindowTitle(title)
        self.setGeometry(top, left, self.width, self.height)
        self.setFixedSize(self.width, self.height)
        
        # Get data class [although this]
        self.inputdata = DATASETTINGS()

        # Run GUI
        self.MyUI()

    def MyUI(self):
        
        # Add plot, entire time course
        self.canvas_tc = PlotCanvas(self,
                                    width=(self.width * .85) / 100,
                                    height=(self.height * .5) / 100)
        self.canvas_tc.move(self.width / 8, 0)
        self.addToolBar(QtCore.Qt.BottomToolBarArea,
                        NavigationToolbar(self.canvas_tc, self))

        # Add plot, trial by trial
        self.canvas_tr = PlotCanvasTrials(self,
                                          width=(self.width * .85) / 100,
                                          height=(self.height * .3) / 100)
        self.canvas_tr.move(self.width / 8, (self.height / 2) + self.height / 10)

        # Add buttons > time course [tc plot]
        x_left = self.width / 60
        y_start = self.width / 60
        self.makebutton("Get data", x=x_left, y=y_start,
                        do_action=self.get_data)
        self.makebutton("Filter", x=x_left, y=y_start * 2,
                        do_action=self.canvas_tc.dofilter)
        self.makebutton("TR corr", x=x_left, y=y_start * 3,
                        do_action=0)
       
        #--peaks
        self.nameLabel = QLabel('Peak frequency:', self)
        self.nameLabel.move(x_left,y_start * 5)
        
        self.makebutton("Plot peaks", x=x_left * 7, y=y_start * 5,
                        do_action=self.canvas_tc.get_peaks)

        mySlider = QSlider(Qt.Horizontal, self)
        mySlider.setGeometry(x_left, y_start * 6, 120, 30)
        mySlider.setMinimum(4)
        mySlider.setMaximum(8)
        mySlider.setValue(6)
        mySlider.setTickInterval(1)
        mySlider.setTickPosition(QSlider.TicksBelow)
        
        
        mySlider.valueChanged[int].connect(self.canvas_tc.slidervalue)
        mySlider.setSingleStep(10)
        #---
        
        self.makebutton("Insert", x=x_left, y=y_start * 8,
                        do_action=self.insertpeak)
        self.makebutton("Remove", x=x_left, y=y_start * 9,
                        do_action=self.removepeaks)
        self.makebutton("Interpolate", x=x_left, y=y_start * 10,
                        do_action=self.canvas_tc.interpol)
        self.makebutton("Save", x=x_left, y=y_start * 12,
                        do_action=self.canvas_tc.save)

        # Add buttons > conversion
        #self.makebutton("Convert EEG", x=x_left, y=y_start * 12,
        #                do_action=self.readineeg)
        #self.makebutton("Convert ACQ", x=x_left, y=y_start * 13,
        #                do_action=self.readinacq)

        # Add buttons > trials [tr canvas]
        self.makebutton("Get trials", x=x_left, y=((self.height / 2) + self.height / 10) + y_start,
                        do_action=self.get_trials)  # self.canvas_tr.get_trials)
        self.makebutton(">", x=x_left, y=((self.height / 2) + self.height / 10) + y_start * 2,
                        do_action=self.canvas_tr.trialup)
        self.makebutton("<", x=x_left, y=((self.height / 2) + self.height / 10) + y_start * 3,
                        do_action=self.canvas_tr.trialdown)
        self.makebutton("Accept", x=x_left, y=((self.height / 2) + self.height / 10) + y_start * 4,
                        do_action=self.canvas_tr.accepttrial)
        self.makebutton("Reject", x=x_left, y=((self.height / 2) + self.height / 10) + y_start * 5,
                        do_action=self.canvas_tr.rejecttrial)
        self.makebutton("Save", x=x_left, y=((self.height / 2) + self.height / 10) + y_start * 6,
                        do_action=self.canvas_tr.save)

        self.fileLabel = QLabel('file name', self)
        self.fileLabel.move(x_left,((self.height / 2) + self.height / 10) + y_start * 8)

    def change_label(self):

        nameoffile = self.canvas_tc.filename[0].rstrip(os.sep)
        head, tail = os.path.split(nameoffile)
        self.fileLabel.setText(tail)

    def makebutton(self, tekst, x, y, do_action):
        """
        Function to make a button to not repeate code
        """
        button = QPushButton(tekst, self)
        button.move(x, y)
        if do_action:
            button.clicked.connect(do_action)
        if tekst == "Reject":
            button.setStyleSheet("color:rgb(125,60,60)")
        elif tekst == "Accept":
            button.setStyleSheet("color:rgb(60,125,60)")

    # Run trial analysis in other class, this contructions allows the
    # two PlotCanvas' to communicate [tc = time course, tr =  trial]
    def get_data(self):

        # Clear the plots in case of reopening a new file
        self.canvas_tc.axes.cla()
        self.canvas_tc.draw()
        self.canvas_tr.axes.cla()
        self.canvas_tr.draw()

        # Get data
        try:
            self.canvas_tc.get_data()
            self.change_label()
        except:
            "No file was selected"
            

    def get_trials(self):
        try:
            self.canvas_tc.peaks
            self.canvas_tr.get_trials()
            self.canvas_tr.update_trials(self.canvas_tc.dat, self.canvas_tc.peaks, self.canvas_tc.ibichannel)
        except AttributeError:
            print("Run peak detection first, please.")

    def insertpeak(self):
        try:
            self.canvas_tc.insertpeak()
            self.canvas_tr.update_trials(self.canvas_tc.dat, self.canvas_tc.peaks, self.canvas_tc.ibichannel)
        except:
            pass

    def removepeaks(self):
        try:
            self.canvas_tc.removepeaks()
            self.canvas_tr.update_trials(self.canvas_tc.dat, self.canvas_tc.peaks, self.canvas_tc.ibichannel)
        except:
            pass

    # Convert EEG files
    def readineeg(self):

        # Pop up menu asking for the file
        filename = QFileDialog.getOpenFileName(self,
                                               'Open a data file', '.', 'vhdr files (*.vhdr);;All Files (*.*)')

        try:
            filename[0]
            raw = mne.io.read_raw_brainvision(filename[0])
            raw.load_data()  # preload so you can resample
            raw.resample(sfreq=self.inputdata.Hz)  # resample to 50 HZ
            dat = raw.get_data()  # extract data

            # Save each channel in a txt file
            for c, val in enumerate(raw.ch_names):
                thefile = open(filename[0][:-5] + "_" + val + ".txt", 'w')
                for item in dat[c]:
                    thefile.write("%s\n" % item)
                thefile.close()

                # Get stimulus events
            events = mne.events_from_annotations(raw)
            thefile = open(filename[0][:-5] + "_stimcodes.txt", 'w')
            for c, val in enumerate(list(events[0])):
                if val[0] > 0:
                    thefile.write("%s\n" % val[0])
            thefile.close()
            print("files coverted!")

        except:

            pass

    # convert ACQ files
    def readinacq(self):

        # Pop up menu asking for the file
        filename = QFileDialog.getOpenFileName(self,
                                               'Open a data file', '.', 'acq files (*.acq);;All Files (*.*)')
        try:
            raw = bioread.read_file(filename[0])

            # print channel names
            for c, val in enumerate(list(raw.named_channels.keys())):
                print("Name of channel " + str(c + 1) + " is: " + val)
            # get name of the event channel [user specific]
            text, okPressed = QInputDialog.getText(self, "Get text", "Name of event channel (see printed options):",
                                                   QLineEdit.Normal, "")
            if okPressed and text != '':
                eventname = text
                #
            channel_names = list(raw.named_channels.keys())
            for c, val in enumerate(channel_names):
                thefile = open(filename[0][:-5] + "_" + val + ".txt", 'w')
                res_dat = raw.channels[c].data[::int(raw.samples_per_second / self.inputdata.Hz)]  # resample to 50 HZ
                for item in res_dat:
                    thefile.write("%s\n" % item)
                thefile.close()

                # Get stimulus events
            eventchannel = [c for c, val in enumerate(channel_names) if val == eventname]
            stimcodes = raw.channels[eventchannel[0]].data
            stimcodes = [c for c, val in enumerate(np.diff(stimcodes)) if val == 5]
            stimcodes = [int(val / int(raw.samples_per_second / self.inputdata.Hz)) for val in stimcodes]
            thefile = open(filename[0][:-5] + "_stimcodes.txt", 'w')
            for item in stimcodes:
                thefile.write("%s\n" % item)
            thefile.close()

            print("files coverted :)")

        except:

            pass


# ------------------------------------------------------------------------------

# Canvas for plotting
class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=5, dpi=100):

        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.set_facecolor([.92, .92, .92])  # Color of the plot canvas
        self.axes = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        # Get data classes
        self.DA = DataAnalysis()
        self.inputdata = DATASETTINGS()
        
        self.peakdis=60 # default

    def slidervalue(self, value):
        
        self.peakdis=value*10
    
    def get_data(self):

        # Get datafile
        self.filename = QFileDialog.getOpenFileName(self,
                                                    'Open a data file', '.', 'TXT files (*.txt);;All Files (*.*)')
        
        
        # Get HR time course
        self.dat = self.openfile(self.filename[0])
        self.dat = self.dat - np.mean(self.dat)
        self.npdat = np.array(self.dat)  # For peak detection works only on array's...

        # In case a peak file already excists, open it and plot it
        try:

            # Open peak file
            path, file = os.path.split(self.filename[0])
            peakfile = os.path.join(path + "/PulseCor_output/" + file[:-4] + "_peaks.txt")
            self.peaks = self.openfile(peakfile)

            # Get IBI
            self.ibichannel = self.DA.get_ibi(self.dat, self.peaks)

            # Plot peaks and IBI
            self.plotall()

        except:

            # Plot data only for first time
            try:
                del self.peaks
                del self.ibichannel
            except:
                pass    
            
            self.axes.cla()
            self.axes.plot(self.dat, c='C0')
            self.draw()

    def dofilter(self):

        # Check if there is data
        try:
            self.dat
        
        except:
            
            print("Please load data first")
        
        else:
            
            # Apply filter
            lowcut, ok = QInputDialog.getInt(self, "Input", "Lowcut:", 1, 0, 20000, 1)
            
            if ok:
                
                highcut, ok = QInputDialog.getInt(self, "Input", "Highcut:", 10, 0, 20000, 1)

                if ok:
                
                    self.dat = self.DA.butter_bandpass_filter(self.dat, lowcut, highcut, self.inputdata.Hz, order=5)
        
                    # Check if peaks was already loaded from before
                    try:
                        self.peaks
                
                    except:
                        
                        # Clear plot
                        self.axes.cla()
        
                        # Plot for first time
                        self.axes.plot(self.dat, c='C0')
                        self.draw()
                    
                    else:
                        
                        # Replot
                        self.remove()
                        self.plotall()


    def get_peaks(self):

        # Get input
        #peakdis, ok = QInputDialog.getInt(self, "Input", "Peak distance (40-90):", 60, 0, 20000, 1)

        
    
        # Get peaks
        self.peaks = self.DA.runpeaks(self.dat, self.peakdis)

        # Get IBI
        self.ibichannel = self.DA.get_ibi(self.dat, self.peaks)

        # Plot
        self.remove() # in case of replotting
        self.plotall()
    
    
    def insertpeak(self):

        try:

            # get xaxis
            x = self.axes.get_xlim()
            x = list(x)

            # get data in the zoomed window
            zoomdat = self.dat[int(x[0]):int(x[1])]

            # get max value within window
            maxpulse = max(zoomdat)
            temppeak = [c for c, val in enumerate(zoomdat) if val == maxpulse]

            # check if peak is already there [+/- 5]
            peakthere = [val for val in self.peaks if
                         val > int(x[0] + temppeak[0]) - 5 and val < int(x[0] + temppeak[0]) + 5]

            # in case it is not add
            if not peakthere:
                # Add peak
                self.peaks.append(int(x[0] + temppeak[0]))

                # Make sure peaks below 0 are removed
                self.peaks = [val for val in self.peaks if val > 0]

                # Recalc IBI
                self.peaks.sort()
                self.ibichannel = self.DA.get_ibi(self.dat, self.peaks)

                # Replot
                self.remove()
                self.plotall()

        except:
            pass

    def removepeaks(self):

        try:

            # get xaxis
            x = self.axes.get_xlim()
            x = list(x)

            # find the peak winthin zoomed in window
            temppeak = [val for c, val in enumerate(self.peaks) if val < x[1] and val > x[0]]

            # if there is one
            if temppeak:
                # get the exact peak in the list
                peakthere = [val for val in self.peaks if val > temppeak[0] - 5 and val < temppeak[0] + 5]

                # remove from peak list
                self.peaks.remove(peakthere[0])

                # Recalc IBI
                self.peaks.sort()
                self.ibichannel = self.DA.get_ibi(self.dat, self.peaks)

                # Replot
                self.remove()
                self.plotall()

        except:
            pass

    def interpol(self):

        # get xaxis
        x = self.axes.get_xlim()
        x = list(x)

        # collect area
        baddata = [int(x[0]), int(x[1])]

        # only interpolate when data is within the window to prevent errors when
        # zoomed in on the edges
        if x[0] > 0 and x[1] < len(self.dat):

            # find the peaks winthin zoomed in window
            temppeak = [val for c, val in enumerate(self.peaks) if val < x[1] and val > x[0]]

            # loop over peaks to remove them
            for p in temppeak:
                # get the exact peak in the list
                peakthere = [val for val in self.peaks if val > p - 5 and val < p + 5]

                # remove from peak list
                self.peaks.remove(peakthere[0])

            # find the closest peak before and after to be interpolated part
            takeClosest = lambda num, collection: min(collection, key=lambda x: abs(x - num))
            st = int(takeClosest(baddata[0], self.peaks))
            en = int(takeClosest(baddata[1], self.peaks))

            # index them
            i_st = [c for c, val in enumerate(self.peaks) if val == st]
            i_en = [c for c, val in enumerate(self.peaks) if val == en]

            # get window before and after [3 sec window?]
            win_st = int(i_st[0] - (self.inputdata.Hz * 3))
            win_en = int(i_en[0] + (self.inputdata.Hz * 3))
            
            
            
            # only when there is a window
            #if win_st > 0 and win_en < len(self.dat):
            if win_st < 0:
                win_st=0
            if win_en > len(self.peaks):
                win_en = len(self.peaks)
                
            print(win_en)
            print(i_en[0])
            print(win_en-i_en[0])
            
            # Enough space on both sides
            if (i_st[0]-win_st) > 0 and (win_en-i_en[0]) > 1:
                
                # calculate the average bmp [thus step size for intyerpolation]
                intval_st = int(np.mean(np.diff(self.peaks[win_st:i_st[0]])))
                intval_en = int(np.mean(np.diff(self.peaks[i_en[0]:win_en])))
                intval = (intval_st + intval_en) / 2
           
            # Not enough space in the beginning, use the end
            elif (i_st[0]-win_st) == 0:
                
                # calculate the average bmp [thus step size for intyerpolation]
                intval = int(np.mean(np.diff(self.peaks[i_en[0]:win_en])))
                
            # Not enough space at the end, use the beginning
            elif (win_en-i_en[0]) == 1:
                
                # Calculate the average bmp [thus step size for interpolation]
                intval = int(np.mean(np.diff(self.peaks[win_st:i_st[0]])))
    
            # calculate how many peaks to insert?
            npeaks = int(np.floor((en - st) / intval)) + 1

            # interpolate
            newpeaks = [int(val) for val in np.arange(st, en, (en - st) / npeaks)]
            self.peaks.extend(newpeaks)
            self.peaks.sort()

            # Recalc IBI
            self.ibichannel = self.DA.get_ibi(self.dat, self.peaks)

            # Replot
            self.remove()
            self.plotall()

    def remove(self):
        for line in self.axes.lines:
            line.remove()

    def plotall(self):

        # Plot [HR, peaks, IBI]
        self.axes.plot(self.peaks, [self.npdat[c] for c, num in enumerate(self.peaks)], '|', color=(0.9, 0.9, 0.9),
                       markersize=500)
        self.axes.plot(self.dat, c='C0')
        self.axes.plot(self.ibichannel, 'y')
        self.draw()

    def openfile(self, file):

        # Extract data from the file
        dat = []
        with open(file) as f:
            for line in f:
                line = line.split('\n')
                dat.append(float(line[0]))
        return dat

    def save(self):

        try:

            # Split path and file for saving
            path, file = os.path.split(self.filename[0])

            # Create output directory
            if not os.path.exists(os.path.join(path, 'PulseCor_output')):
                os.makedirs(os.path.join(path, 'PulseCor_output'))

            # Sort [because new peaks are added at the end of the list]
            self.peaks.sort()

            # Create file and save in PulseCor_output directory

            newname = os.path.join(path, "PulseCor_output", file[:-4] + "_peaks.txt")
            thefile = open(newname, 'w')
            for item in self.peaks:
                thefile.write("%s\n" % item)
            thefile.close()
            
            newname = os.path.join(path, "PulseCor_output", file[:-4] + "_IBI.txt")
            thefile = open(newname, 'w')
            for item in self.ibichannel:
                thefile.write("%s\n" % item)
            thefile.close()

        except:
            pass


# Canvas for plotting
class PlotCanvasTrials(FigureCanvas):

    def __init__(self, parent=None, width=5, height=5, dpi=100):

        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.set_facecolor([.92, .92, .92])  # Color of the plot canvas
        self.axes = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        # Extra settings
        self.trnum = 0
        self.trialsaccepted = []

        # Get data classes
        self.DA = DataAnalysis()
        self.inputdata = DATASETTINGS()

    def get_trials(self):

        # Get stimulus onset times
        self.stimonset = self.openfile()

        # Get baseline and stimulus duration in seconds
        self.inputdata.baseline, ok = QInputDialog.getInt(self, "Input", "Baseline in sec:", self.inputdata.baseline, 0,
                                                          20000, 1)
        self.inputdata.stimduur, ok = QInputDialog.getInt(self, "Input", "Stim duration in sec:",
                                                          self.inputdata.stimduur, 0, 20000, 1)

        # Make vector with 1's for accepted trials
        self.trialsaccepted = np.ones(len(self.stimonset))

    def update_trials(self, dat, peaks, ibi):
        """
        Fuction is neccessary to make the 2 plots communicate
        """
        self.dat = dat
        self.peaks = peaks
        self.ibichannel = ibi

        self.plot_trials()

    def plot_trials(self):

        # Get first trial
        tr_onset = int(self.stimonset[self.trnum] - self.inputdata.Hz * self.inputdata.baseline)
        adj_stimduur = self.inputdata.Hz * self.inputdata.stimduur
        tr_offset = int(self.stimonset[self.trnum] + adj_stimduur)

        # Cut out data
        self.trialdat = self.dat[tr_onset:tr_offset]

        # Peaks
        self.trialpeaks = [val for val in self.peaks if val > tr_onset]
        self.trialpeaks = [val for val in self.trialpeaks if val < tr_offset]
        self.trialpeaks = [val - tr_onset for val in self.trialpeaks]

        # IBI
        self.trialibi = self.ibichannel[tr_onset:tr_offset]

        # Clear plot
        self.axes.cla()

        # Plot first trial
        if self.trialsaccepted[self.trnum] == 1:
            self.axes.plot(self.trialdat, c='C0')
        else:
            self.axes.plot(self.trialdat, c='gray')
        
        self.axes.plot(self.trialpeaks, [self.trialdat[c] for c, num in enumerate(self.trialpeaks)], '|',
                       color=(0.9, 0.9, 0.9), markersize=500)
        if self.trialsaccepted[self.trnum] == 1:
            self.axes.plot(self.trialibi, 'y')
        else:
            self.axes.plot(self.trialibi, 'gray')
        self.draw()

    def trialup(self):

        try:

            self.trnum += 1

            if self.trnum >= len(self.stimonset):
                self.trnum = len(self.stimonset) - 1

            self.plot_trials()

        except AttributeError:

            pass

    def trialdown(self):

        try:

            self.trnum -= 1

            if self.trnum < 0:
                self.trnum = 0

            self.plot_trials()

        except AttributeError:

            pass

    def accepttrial(self):

        try:

            self.trialsaccepted[self.trnum] = 1
            self.plot_trials()

        except IndexError:

            pass

    def rejecttrial(self):

        try:

            self.trialsaccepted[self.trnum] = 0
            self.plot_trials()

        except IndexError:

            pass

    def openfile(self):

        # Get datafile
        self.filename = QFileDialog.getOpenFileName(self,
                                                    'Open a data file', '.', 'TXT files (*.txt);;All Files (*.*)')

        try:
            # Extract data from the file
            dat = []
            with open(self.filename[0]) as f:
                for line in f:
                    line = line.split('\n')
                    dat.append(float(line[0]))
            return dat

        except:

            pass

    def save(self):

        try:

            path, file = os.path.split(self.filename[0])

            # Create output directory
            if not os.path.exists(os.path.join(path, 'PulseCor_output')):
                os.makedirs(os.path.join(path, 'PulseCor_output'))

            # Create file and save in PulseCor_output directory
            newname = os.path.join(path, "PulseCor_output", file[:-4] + "_acceptedtrials.txt")
            thefile = open(newname, 'w')
            for item in self.trialsaccepted:
                thefile.write("%s\n" % int(item))
            thefile.close()

        except:

            pass


# ------------------------------------------------------------------------------

class DataAnalysis():

    def __init__(self):

        self.inputdata = DATASETTINGS()

    def runpeaks(self, dat, peakdis):

        """
        https://dsp.stackexchange.com/questions/20231/good-way-to-detect-pulse-with-known-width-with-background-noise/21839#21839
        """

        # works only on array data
        data = np.array(dat)

        # get min max depending on filter window
        max_data = maximum_filter(data, peakdis)
        min_data = -maximum_filter(-data, peakdis)

        # select places where we detect maximum but not minimum -> we dont want long plateaus
        peak_mask = np.logical_and(max_data == data, min_data != data)

        # select peaks where we have enough elevation [REMOVED FOR ECG, STILL WORKS FOR PULSE]
        # threshold=1 # this can be changed, but for now always look for 1 peak in a given window
        # peak_mask = np.logical_and(peak_mask, max_data - min_data > threshold)

        # a trick to convert True to 1, False to -1
        peak_mask = peak_mask * 2 - 1

        # select only the up edges to eliminate multiple maximas in a single peak
        peak_mask = np.correlate(peak_mask, [-1, 1], mode='same') == 2
        peaks = np.where(peak_mask)[0]
        peaks = peaks.tolist()  # make list

        return peaks

    def get_ibi(self, dat, peaks):

        peaks = [int(val) for val in peaks]

        ibi = np.diff(peaks)
        ibichannel = np.zeros(len(dat))
        for c, val in enumerate(peaks):

            # Add the beginning
            if c == 0:
                ibichannel[0:peaks[c]] = ibi[c]

            # Add IBI inbetween peaks
            if c + 1 < len(peaks):
                ibichannel[peaks[c]:peaks[c + 1]] = ibi[c]

        # Add the end
        ibichannel[peaks[c]:] = ibi[-1]

        # Adjust range
        ibichannel = np.interp(ibichannel, (ibichannel.min(), ibichannel.max()), (min(dat), max(dat)))

        # zero center
        ibichannel = ibichannel - np.mean(ibichannel)

        return ibichannel

    def butter_bandpass(self, lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a

    def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=5):
        b, a = self.butter_bandpass(lowcut, highcut, fs, order=order)
        y = filtfilt(b, a, data, method='pad')
        return y


# RUN
# ------------------------------------------------------------------------------

app = QApplication(sys.argv)
window = Window()
window.show()
app.exec()

