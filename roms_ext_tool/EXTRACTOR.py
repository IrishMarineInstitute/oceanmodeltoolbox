"""
! ********************************************************************
! *                                                                  *
! *         Irish Marine Institute ROMS Data Extractor               *
! *                         										 *
! * Developer: Diego Pereiro          								 *
! * Release  : v1.1                      							 *
! * Date     : January 2020                      					 *
! *                         										 *
! ********************************************************************
!
!  Created by:            Diego Pereiro
!  Created on:            15 January 2020
!  Last Modified on:      11 February 2020
!
!
!
! **********************************************************************
! **********************************************************************
! **                      Copyright (c) 2016                          **
! **               The Marine Institute, Ireland                      **
! **********************************************************************
! **                                                                  **
! ** This Software is open-source and licensed under the following    **
! ** conditions as stated by MIT/X License:                           **
! **                                                                  **
! **  (See http://www.opensource.org/licenses/mit-license.php ).      **
! **                                                                  **
! ** Permission is hereby granted, free of charge, to any person      **
! ** obtaining a copy of this Software and associated documentation   **
! ** files (the "Software"), to deal in the Software without          **
! ** restriction, including without limitation the rights to use,     **
! ** copy, modify, merge, publish, distribute, sublicense,            **
! ** and/or sell copies of the Software, and to permit persons        **
! ** to whom the Software is furnished to do so, subject to the       **
! ** following conditions:                                            **
! **                                                                  **
! ** The above copyright notice and this permission notice shall      **
! ** be included in all copies or substantial portions of the         **
! ** Software.                                                        **
! **                                                                  **
! ** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  **
! ** EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE           **
! ** WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE  **
! ** AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  **
! ** HOLDERS BE LIABLE FOR ANY CLAIMS, DAMAGES OR OTHER LIABILITIES,  **
! ** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     **
! ** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR    **
! ** OTHER DEALINGS IN THE SOFTWARE.                                  **
! **                                                                  **
! ** The most current official versions of this Software and          **
! ** associated tools and documentation are available from the        ** 
! ** authors by e-mail:                                               **
! **                                                                  **
! **        diego.pereiro@marine.ie                                   **
! **                                                                  **
! ** We ask that users make appropriate acknowledgement of            **
! ** The Marine Institute, Ireland,                                   **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies.                                        **
! **********************************************************************
! **********************************************************************
"""

import warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import netCDF4
import scipy
import cftime
import re
from tkinter import ttk
from datetime import datetime, timedelta
from math import ceil, floor

class Root(tk.Tk):  
    
    def __init__(self):        
        super(Root, self).__init__()
        # Window title
        self.title("Marine Institute ROMS Data Extractor")
        # Window size                
        self.minsize(50, 50); self.resizable(0, 0)
                      
        self.numtimes = 0
        self.choices = ()
        
        """ Read config file """
        self.readconfig() 
        self.step_config = max(self.tstep_config, self.mstep_config)
        
        """ Read coastline """
        self.coast_x = []
        self.coast_y = []
        infile = open("mundo_h.dat", "r")
        for line in infile:
            lst = line.split()
            self.coast_x.append(float(lst[0]))
            self.coast_y.append(float(lst[1]))  
            
        """ Returning variables """
        # Start times
        self.t0 = []
        # End times
        self.t1 = []        
        # User-selected variables
        self.userkeys = []
        # NetCDF variable descriptions
        self.descriptions = []
        # Minimum of colorbar range
        self.mincolran = []
        # Maximum of colorbar range
        self.maxcolran = []
        # Color Maps
        self.colors = []        
        # Color Scales
        self.scales = []
        
        """ Variable dictionary """
        self.numvars = 0
        self.vardict = self.varStruct(self.choices)
        self.keylist = list(self.vardict)
        self.long = [item[2] for item in self.vardict.values()]
        self.cscale = [item[7] for item in self.vardict.values()]
        self.caxis = [item[8] for item in self.vardict.values()]
        self.cmap = [item[9] for item in self.vardict.values()]
        
        """ Browse frame """
        # Frame
        self.browseFrame = tk.LabelFrame(self, text="Choose a file", padx=5, pady=5)
        self.browseFrame.grid(row=0, column=0, sticky=tk.W + tk.E, padx=10, pady=10)        
        # Browse button
        self.browseButton = tk.Button(self.browseFrame, text="Select output directory", command=self.fileDialog)
        self.browseButton.grid(row=0, column=0, padx=5, pady=5)
        # Browse entry
        self.browseEntry = tk.Entry(self.browseFrame, width=50)
        self.browseEntry.grid(row=0, column=1, padx=5, pady=5) 
        # Filename label
        self.fileLabel = tk.Label(self.browseFrame, text="Output file name")
        self.fileLabel.grid(row=0, column=2, padx=5, pady=5)
        # Filename entry
        self.fileEntry = tk.Entry(self.browseFrame, width=20)
        self.fileEntry.grid(row=0, column=3, padx=5, pady=5)
        
        """ Mode frame """
        # Frame
        self.modeFrame = tk.LabelFrame(self, text="Select mode", padx=5, pady=5)
        self.modeFrame.grid(row=0, column=1, sticky=tk.W + tk.E, padx=10, pady=10)
        # Mode Radiobutton          
        self.mode = tk.IntVar(value=1)
        self.mode1 = tk.Radiobutton(self.modeFrame, text="Time averages", variable=self.mode, value=1)
        self.mode1.grid(row=0, column=0, sticky=tk.W, padx=40)
        self.mode2 = tk.Radiobutton(self.modeFrame, text="Time series", variable=self.mode, value=2)
        self.mode2.grid(row=0, column=1, sticky=tk.W, padx=40)
                 
        """ Z-levels frame """        
        # Frame
        self.zFrame = tk.LabelFrame(self, text="Select z-levels", padx=5, pady=25)
        self.zFrame.grid(row=2, column=0, sticky=tk.W + tk.E, padx=10, pady=10)
        # Z-levels instructions
        self.zLabel = tk.Label(self.zFrame, text="comma-separated depths [meter]")
        self.zLabel.grid(row=0, column=0, padx=0)
        # Z-levels entry
        self.zEntry = tk.Entry(self.zFrame, width=80)        
        self.zEntry.grid(row=0, column=1, padx=5, pady=5)
                       
        """ Process dates """
        self.dates = []; self.mon = []
        self.time0 = self.todatetime(self.idate_config)               
        self.time1 = self.todatetime(self.edate_config)               
        self.time1 += timedelta(seconds=self.step_config)
        while self.time0 < self.time1:
            self.dates.append(self.time0)
            self.mon.append(self.time0.year * 12 + self.time0.month)
            self.time0 += timedelta(seconds=self.step_config)
                
        """ Period frame """
        # Frame
        self.periodFrame = tk.LabelFrame(self, text="Time period", padx=5, pady=5)
        self.periodFrame.grid(row=1, column=0, sticky=tk.W + tk.E, padx=10, pady=10) 
        
        # Quick-selection frame
        self.quickFrame = tk.LabelFrame(self.periodFrame, text="Quick selection")
        self.quickFrame.grid(row=0, column=0, padx=5, pady=5)
        # "Every" selection
        self.everyLabel = tk.Label(self.quickFrame, text="Every")
        self.everyLabel.pack(side=tk.LEFT, padx=8, pady=5)
        # Step CheckButton
        self.steplySwitch = tk.IntVar(value=0)
        self.steplyCheckButton = tk.Checkbutton(self.quickFrame, text = "step", variable = self.steplySwitch)
        self.steplyCheckButton.pack(side=tk.LEFT, padx=8, pady=5)        
        # Daily Checkbutton
        self.dailySwitch = tk.IntVar(value=0)
        self.dailyCheckButton = tk.Checkbutton(self.quickFrame, text = "day", variable = self.dailySwitch)
        self.dailyCheckButton.pack(side=tk.LEFT, padx=8, pady=5)
        # Weekly CheckButton
        self.weeklySwitch = tk.IntVar(value=0)
        self.weeklyCheckButton = tk.Checkbutton(self.quickFrame, text = "week", variable = self.weeklySwitch)
        self.weeklyCheckButton.pack(side=tk.LEFT, padx=8, pady=5)
        # Monthly CheckButton
        self.monthlySwitch = tk.IntVar(value=0)
        self.monthlyCheckButton = tk.Checkbutton(self.quickFrame, text = "month", variable = self.monthlySwitch)
        self.monthlyCheckButton.pack(side=tk.LEFT, padx=8, pady=5)
        # "From" label
        self.fromLabel = tk.Label(self.quickFrame, text="from")
        self.fromLabel.pack(side=tk.LEFT, padx=8, pady=5)
        # Initial month drop-down menu
        self.monthsList = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]        
        self.imonth = ttk.Combobox(self.quickFrame, values=self.monthsList, state="readonly", width=5)
        self.imonth.pack(side=tk.LEFT, padx=7, pady=5)
        # Initial year drop-down menu
        self.iyear = ttk.Combobox(self.quickFrame, values=list(range(self.idate_config[0], 1 + self.edate_config[0])), state="readonly", width=5)
        self.iyear.pack(side=tk.LEFT, padx=7, pady=5)
        # "To" label
        self.toLabel = tk.Label(self.quickFrame, text="to")
        self.toLabel.pack(side=tk.LEFT, padx=8, pady=5)
        # End month drop-down menu
        self.emonth = ttk.Combobox(self.quickFrame, values=self.monthsList, state="readonly", width=5)
        self.emonth.pack(side=tk.LEFT, padx=7, pady=5)
        # End year drop-down menu
        self.eyear = ttk.Combobox(self.quickFrame, values=list(range(self.idate_config[0], 1 + self.edate_config[0])), state="readonly", width=5)
        self.eyear.pack(side=tk.LEFT, padx=7, pady=5)        
        
        # Manual frame
        self.manualFrame = tk.LabelFrame(self.periodFrame, text="Manual selection")
        self.manualFrame.grid(row=1, column=0, padx=5, pady=5)
        # "start" label
        self.startLabel = tk.Label(self.manualFrame, text="From")
        self.startLabel.grid(row=0, column=0, padx=5, pady=5)
         # Start drop-down menu
        self.idate = ttk.Combobox(self.manualFrame, values=self.dates, state="readonly", width=37)
        self.idate.grid(row=0, column=1, padx=5)
        # "end" label
        self.endLabel = tk.Label(self.manualFrame, text="to")
        self.endLabel.grid(row=0, column=2, padx=5)
        # End time drop-down menu
        self.edate = ttk.Combobox(self.manualFrame, values=self.dates, state="readonly", width=37)
        self.edate.grid(row=0, column=3, padx=5)        
        # Add button
        self.addPeriodButton = tk.Button(self.manualFrame, text="Add", command=self.addPeriod)
        self.addPeriodButton.grid(row=0, column=6, padx=5)
        # Del button
        self.delTimeButton = tk.Button(self.manualFrame, text="Delete", command=self.delTime)
        self.delTimeButton.grid(row=0, column=7, padx=5)
        # Time period entry field
        self.timeEntry = tk.Entry(self.manualFrame, width=90, state="disabled", \
                                  disabledbackground="white", disabledforeground="black")
        self.timeEntry.grid(row=1, column=0, columnspan=8, sticky=tk.E + tk.W, padx=5, pady=5)
        # Scrollbar
        self.tscrollbar = ttk.Scrollbar(self.manualFrame, command=self.timeEntry.xview, orient=tk.HORIZONTAL)
        self.tscrollbar.grid(row=2, column=0, columnspan=8, sticky=tk.E + tk.W, padx=5)
        self.timeEntry["xscrollcommand"] = self.tscrollbar.set
     
        """ Get input file name matching initial date """
        self.rep = {"%Y": self.idate_config_str[0], \
               "%y": self.idate_config_str[0][2:], \
               "%m": self.idate_config_str[1], \
               "%b": self.monthsList[self.idate_config[1]-1], \
               "%d": self.idate_config_str[2], \
               "%o": str(datetime(self.idate_config[0], self.idate_config[1], self.idate_config[2]).timetuple().tm_yday).zfill(3), \
               "%H": self.idate_config_str[3], \
               "%M": self.idate_config_str[4], \
               "%S": self.idate_config_str[5], \
               "%rt": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.tstep_config, \
                                             self.todatetime(self.idate_config)), \
               "%rm": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.mstep_config, \
                                             self.todatetime(self.idate_config))}
        self.rep = dict((re.escape(kk), vv) for kk, vv in self.rep.items()) 
        pattern = re.compile("|".join(self.rep.keys()))
        f = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.tracers_convention)
        f = netCDF4.Dataset(f, "r")
        
        """ Read coordinates from input file """
        self.lon = f.variables["lon_rho"][:]; Mp, Lp = self.lon.shape
        self.lat = f.variables["lat_rho"][:]
        
        """ Select area frame """
        # Frame 
        self.areaFrame = tk.LabelFrame(self, text="Select area", padx=5, pady=5)
        self.areaFrame.grid(row=1, column=1, rowspan=2, sticky=tk.W + tk.E + tk.S + tk.N, padx=10, pady=10)
        
        # Axis limits    
        self.minx = self.lon.min()
        self.maxx = self.lon.max()
        self.miny = self.lat.min()
        self.maxy = self.lat.max()
        
        # Subset coast
        self.subsetcoast(self.coast_x, self.coast_y)     
                
        # West string
        self.defaultWest = tk.StringVar(value="%d" % 0)
        self.defaultWest.trace_add("write", self.drawArea)       
        # East string        
        self.defaultEast = tk.StringVar(value="%d" % (Lp - 1))
        self.defaultEast.trace_add("write", self.drawArea)
        # South string
        self.defaultSouth = tk.StringVar(value="%d" % 0)
        self.defaultSouth.trace_add("write", self.drawArea)
        # North string
        self.defaultNorth = tk.StringVar(value="%d" % (Mp - 1))
        self.defaultNorth.trace_add("write", self.drawArea)
        
        # Plot map
        self.fig = Figure(figsize=(3.9576, 3), facecolor=(.95, .95, .95))
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlim(left=self.minx, right=self.maxx); self.ax.get_xaxis().set_ticks([])
        self.ax.set_ylim(bottom=self.miny, top=self.maxy);   self.ax.get_yaxis().set_ticks([])
        self.drawshore(self.ax, self.coast_x, self.coast_y)           
        
        # North entry
        self.northEntry = tk.Entry(self.areaFrame, width=7, textvariable=self.defaultNorth)
        self.northEntry.grid(row=0, column=1, padx=5, pady=5)
        # West entry
        self.westEntry = tk.Entry(self.areaFrame, width=7, textvariable=self.defaultWest)
        self.westEntry.grid(row=1, column=0, padx=5, pady=5)              
        # East entry
        self.eastEntry = tk.Entry(self.areaFrame, width=7, textvariable=self.defaultEast)
        self.eastEntry.grid(row=1, column=2, padx=5, pady=5)
        # South entry
        self.southEntry = tk.Entry(self.areaFrame, width=7, textvariable=self.defaultSouth)
        self.southEntry.grid(row=2, column=1, padx=5, pady=5)      
         
        # Set canvas
        self.canvas = FigureCanvasTkAgg(self.fig, self.areaFrame)        
        self.canvas.get_tk_widget().grid(row=1, column=1)
        
        self.line = []          
                
        """ Variable frame """
        # Frame
        self.varFrame = tk.LabelFrame(self, text="Select variables", padx=5, pady=5)
        self.varFrame.grid(row=3, column=0, columnspan=2, sticky=tk.W + tk.E, padx=10, pady=10)    
        # "name" label
        self.nameLabel = tk.Label(self.varFrame, text="name")
        self.nameLabel.grid(row=0, column=0, padx=0)
        # Key drop-down menu
        self.keystring = tk.StringVar(value="")
        self.key = ttk.Combobox(self.varFrame, textvariable=self.keystring, \
                                values=self.keylist, postcommand=self.rmVar, state="readonly", width=30)
        self.key.bind("<<ComboboxSelected>>", self.TextBoxUpdate)
        self.key.grid(row=0, column=1, padx=5)
        # "title" label
        self.titleLabel = tk.Label(self.varFrame, text="title")
        self.titleLabel.grid(row=0, column=2, padx=0)        
        # Title entry field
        self.title = tk.Entry(self.varFrame, width=40)
        self.title.grid(row=0, column=3, padx=5)
        # Minimum color label
        self.minlabel = tk.Label(self.varFrame, text="min")
        self.minlabel.grid(row=0, column=4, padx=0)
        # Minimum color entry
        self.minimo = tk.Entry(self.varFrame, width=5)
        self.minimo.grid(row=0, column=5, padx=5)        
        # Maximum color label
        self.maxlabel = tk.Label(self.varFrame, text="max")
        self.maxlabel.grid(row=0, column=6, padx=0)
        # Minimum color entry
        self.maximo = tk.Entry(self.varFrame, width=5)
        self.maximo.grid(row=0, column=7, padx=5)  
        # Color Map label
        self.maplabel = tk.Label(self.varFrame, text="colormap")
        self.maplabel.grid(row=0, column=8, padx=0)
        # Color Map drop-down menu
        self.colormaps = ttk.Combobox(self.varFrame, values=plt.colormaps(), state="readonly")
        self.colormaps.grid(row=0, column=9, padx=5)
        # Color Scale label
        self.scalelabel = tk.Label(self.varFrame, text="scale")
        self.scalelabel.grid(row=0, column=10, padx=0)
        # Color Scale radio buttons
        self.scale = tk.IntVar()
        self.lin = tk.Radiobutton(self.varFrame, text="linear", variable=self.scale, value=1)
        self.lin.grid(row=0, column=11, padx=5)
        self.log = tk.Radiobutton(self.varFrame, text="log", variable=self.scale, value=2)
        self.log.grid(row=0, column=12, padx=5)
        # Variable entry field
        self.varEntry = tk.Entry(self.varFrame, width=160, state="disabled", \
                                 disabledbackground="white", disabledforeground="black")
        self.varEntry.grid(row=1, column=0, columnspan=15, sticky=tk.E + tk.W, padx=5, pady=5)
        # Scrollbar
        self.vscrollbar = ttk.Scrollbar(self.varFrame, command=self.varEntry.xview, orient=tk.HORIZONTAL)
        self.vscrollbar.grid(row=2, column=0, columnspan=15, sticky=tk.E + tk.W, padx=5)
        self.varEntry["xscrollcommand"] = self.vscrollbar.set  
        # Add button
        self.addVarButton = tk.Button(self.varFrame, text="Add", command=self.addVar)
        self.addVarButton.grid(row=0, column=13, padx=5) 
        # Del button
        self.delVarButton = tk.Button(self.varFrame, text="Delete", command=self.delVar)
        self.delVarButton.grid(row=0, column=14, padx=5)
        
        """ Plotting frame """
        # Frame
        self.plotFrame = tk.Frame(self, padx=5, pady=5)
        self.plotFrame.grid(row=4, column=0, sticky=tk.W, padx=0)
        # MakePlot check button
        self.makeplot = tk.IntVar(value=0)
        self.mapCheckButton = tk.Checkbutton(self.plotFrame, text = "Create plots", variable = self.makeplot)
        self.mapCheckButton.grid(row=0, column=0, sticky=tk.W, padx=5)        
        # Image format         
        self.fmt = tk.IntVar(value=3)
        self.eps = tk.Radiobutton(self.plotFrame, text="EPS", variable=self.fmt, value=1)
        self.eps.grid(row=0, column=1, sticky=tk.W, padx=5)
        self.pdf = tk.Radiobutton(self.plotFrame, text="PDF", variable=self.fmt, value=2)
        self.pdf.grid(row=0, column=2, sticky=tk.W, padx=5)
        self.png = tk.Radiobutton(self.plotFrame, text="PNG", variable=self.fmt, value=3)
        self.png.grid(row=0, column=3, sticky=tk.W, padx=5)
              
        """ OK frame """
        self.ok = tk.Button(self, text="OK", command=self.accept)   
        self.ok.grid(row=5, column=0, columnspan=2, pady=5)    
        
        """ Display frame """
        # Frame
        self.displayFrame = tk.LabelFrame(self, text="Display", padx=5, pady=5)
        self.displayFrame.grid(row=6, column=0, columnspan=2, sticky=tk.W + tk.E, padx=10, pady=10)
        # Summary display
        self.summaryDisplay = tk.Text(self.displayFrame, height=10, width=140, padx=5, pady=5)
        self.summaryDisplay.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)           
        # Scrollbar
        self.scrollbar = ttk.Scrollbar(self.displayFrame, command=self.summaryDisplay.yview)
        self.scrollbar.grid(row=0, column=1, sticky="nsew")
        self.summaryDisplay["yscrollcommand"] = self.scrollbar.set   
        
        """ Yes/No frame """
        self.usersInput = tk.IntVar()
        # Frame
        self.inputFrame = tk.Frame(self, padx=5, pady=5)
        self.inputFrame.grid(row=7, column=0, columnspan=2, padx=0)
        # Accept button
        self.acceptButton = tk.Button(self.inputFrame, text="Yes", \
                                      foreground="green", command=lambda: self.usersInput.set(1))
        self.acceptButton.grid(row=0, column=0, padx=5, pady=5)
        # Cancel button
        self.cancelButton = tk.Button(self.inputFrame, text="No", \
                                      foreground="red", command=lambda: self.usersInput.set(0))
        self.cancelButton.grid(row=0, column=1, padx=5, pady=5)  
        
        # Upper-left icon
        self.wm_iconbitmap("icon.ico") 
        
    def drawshore(self, ax, x, y, nfill=100):
        """ Draw the coastline """
        import numpy as np
        # Draw the coastline
        ax.plot(x, y, "k-")
        # Search for missing values
        w = np.isnan(x); w = [i for i, val in enumerate(w) if val]
        # Fill  coastal polygons
        for idx in range(min(nfill, len(w)-1)):
            X = x[w[idx]+1:w[idx+1]]
            Y = y[w[idx]+1:w[idx+1]]
            ax.fill(X, Y, ".5")
        
    def subsetcoast(self, x, y):
        import numpy as np
        from math import nan
        # Output initialization
        x_coast = [nan]; y_coast = [nan]
        # Search for missing values
        w = np.isnan(x); w = [i for i, val in enumerate(w) if val]
        for idx in range(len(w)-1):
            # Subset coastal polygon
            X = x[w[idx]+1:w[idx+1]]; Y = y[w[idx]+1:w[idx+1]]
            for i, j in zip(X, Y):
                if (i > self.minx and i < self.maxx and j > self.miny and j < self.maxy):
                    x_coast += X; x_coast.append(nan)
                    y_coast += Y; y_coast.append(nan); break
        self.coast_x = x_coast; self.coast_y = y_coast
        
    def get_running_number(self, idate, step, date):                    
        return str(int(1 + (date - idate).total_seconds() / step)).zfill(4)
    
    def todatetime(self, timelist):
        # Convert datevec to datetime
        [y, m, d, H, M, S] = timelist
        return datetime(y, m, d, H, M, S)
        
    def readconfig(self):        
        infile = open("config.txt", "r")        
        for line in infile:
            if line[0] != "!":
                entry = line.split("=")[1].strip()
                if line[0] == "T":
                    # TRACERS CONVENTION
                    self.tracers_convention = entry
                    self.momentum_convention = entry
                elif line[0] == "M":
                    # MOMENTUM CONVENTION
                    self.momentum_convention = entry
                elif line[0] == "I":
                    # INITIAL DATE
                    self.idate_config_str = re.findall(r"[\w']+", entry)
                    self.idate_config = [int(i) for i in self.idate_config_str]
                elif line[0] == "E":
                    # END DATE
                    self.edate_config_str = re.findall(r"[\w']+", entry)
                    self.edate_config = [int(i) for i in self.edate_config_str]     
                elif line[0] == "S":
                    # TSTEP
                    self.tstep_config = int(entry)  
                    self.mstep_config = int(entry)
                elif line[0] == "Z":
                    # MSTEP
                    self.mstep_config = int(entry)
                elif line[0] == "O":
                    # OFFSET
                    off = [int(i) for i in re.findall(r"[\w']+", entry)]                    
                    self.offset = self.todatetime(off)                    
                elif line[0] == "V":
                    # VARIABLES
                    self.choices = tuple([int(i) for i in entry.split(",")])
        
    def drawArea(self, *args):
        # If any, remove previous line
        if len(self.line):
            self.line.pop(0).remove()
            
        """ Read user-defined boundaries """
        # Western boundary        
        wb = int(self.westEntry.get()) 
        # Eastern boundary
        eb = int(self.eastEntry.get())
        # Southern boundary
        sb = int(self.southEntry.get())
        # Northern boundary
        nb = int(self.northEntry.get()) 
        
        """ Area vertices """
        # Longitude
        xv = [self.lon[sb, wb], self.lon[nb, wb], self.lon[nb, eb], self.lon[sb, eb], self.lon[sb, wb]]
        # Latitude
        yv = [self.lat[sb, wb], self.lat[nb, wb], self.lat[nb, eb], self.lat[sb, eb], self.lat[sb, wb]]
        
        # Plot user's selection
        self.line = self.ax.plot(xv, yv, "r-", linewidth=2)
        # Update canvas       
        self.canvas.draw()        
        
    def delVar(self):
        # Number of variables before "delete"
        n0 = self.numvars
        # Number of variables after "delete"
        n1 = n0 - 1
        if n0:            
            self.numvars = n1
            # Remove last variable
            del self.userkeys[-1]              
            del self.descriptions[-1]
            del self.mincolran[-1]
            del self.maxcolran[-1]
            del self.colors[-1]
            del self.scales[-1]
            # Modify text in entry field
            self.varEntry.configure(state="normal")
            string = self.varEntry.get()
            substr = string[string.find("%d. KEY" % n0):string.find("%d. KEY" % n1)]            
            self.varEntry.delete(0, tk.END)            
            self.varEntry.insert(tk.END, string.replace(substr, ""))            
            self.varEntry.configure(state="disabled")
    
    def delTime(self):
        # Number of time periods before "delete"
        n0 = self.numtimes
        # Number of time periods after "delete"
        n1 = n0 - 1
        if n0:
            self.numtimes = n1
            # Remove last time period
            del self.t0[-1]        
            del self.t1[-1]            
            # Modify text in entry field
            self.timeEntry.configure(state="normal")
            string = self.timeEntry.get()
            substr = string[string.find("%d. " % n0):string.find("%d. " % n1)]    
            self.timeEntry.delete(0, tk.END)
            self.timeEntry.insert(tk.END, string.replace(substr, ""))
            self.timeEntry.configure(state="disabled")
        
    def TextBoxUpdate(self, event):
        self.title.delete(0, tk.END)
        self.title.insert(0, self.long[self.keylist.index(self.key.get())])
        self.minimo.delete(0, tk.END)
        self.minimo.insert(0, self.caxis[self.keylist.index(self.key.get())][0])
        self.maximo.delete(0, tk.END)
        self.maximo.insert(0, self.caxis[self.keylist.index(self.key.get())][1]) 
        if self.cscale[self.keylist.index(self.key.get())] == "linear":
            self.scale.set(1)
        elif self.cscale[self.keylist.index(self.key.get())] == "log":
            self.scale.set(2)
        self.colormaps.set(self.cmap[self.keylist.index(self.key.get())])
        
    def rmVar(self):
        updtlist = [en for en in self.keylist if en not in self.userkeys]
        self.key["values"] = updtlist     
        
    def addVar(self):        
        # Check if "name" field has been selected. Otherwise, do nothing 
        en = self.key.get()
        if en:
            """ Get inputs """
            # Description
            lg = self.title.get()
            # Color min.
            c0 = self.minimo.get() 
            try:
                float(c0)
            except ValueError:      
                self.summaryDisplay.delete("1.0", tk.END)
                self.summaryDisplay.insert(tk.END, "Warning! Not a valid minimum color range. " + \
                                           "Use a dot (.) as decimal separator\n\n")
                return
            # Color max.
            c1 = self.maximo.get()
            try:
                float(c1)
            except ValueError:  
                self.summaryDisplay.delete("1.0", tk.END)
                self.summaryDisplay.insert(tk.END, "Warning! Not a valid maximum color range. " + \
                                           "Use a dot (.) as decimal separator\n\n")
                return
            # Color map
            cm = self.colormaps.get() 
            # Color scale             
            if self.scale.get() == 1:
                self.scales.append("linear")
            elif self.scale.get() == 2:
                self.scales.append("log")
            """ Update number of variables so far """
            self.numvars = self.numvars + 1    
            """ Input into variable field """
            self.varEntry.configure(state="normal")
            self.varEntry.insert(0, str(self.numvars) + ". " + "KEY: " + en + \
                             ", DESCRIPTION: " + lg + ", Y-MIN: " + c0 + ", Y-MAX: " + c1 + \
                             ", COLORMAP: " + cm + ", Y-SCALE: " + self.scales[-1] + "      ")        
            self.varEntry.configure(state="disabled")
            """ Append new selections """
            self.userkeys.append(en)
            self.descriptions.append(lg)
            self.mincolran.append(float(c0))
            self.maxcolran.append(float(c1))
            self.colors.append(cm)   
        """ Clear fields """
        self.key.set("")
        self.title.delete(0, tk.END)
        self.minimo.delete(0, tk.END)
        self.maximo.delete(0, tk.END)
        self.colormaps.set("")
                
    def addPeriod(self):
        time0 = self.idate.get()
        time1 = self.edate.get()        
        self.numtimes = self.numtimes + 1
        self.timeEntry.configure(state="normal")
        if time0 and time1:
            if time0 != time1:                
                self.timeEntry.insert(0, str(self.numtimes) + ". From " + time0 + " to " + time1)                
            else:
                self.timeEntry.insert(0,  str(self.numtimes) + ". Snapshot " + time0 + "      ")     
            self.t0.append(datetime.strptime(time0, "%Y-%m-%d %H:%M:%S"))
            self.t1.append(datetime.strptime(time1, "%Y-%m-%d %H:%M:%S"))           
        self.timeEntry.configure(state="disabled")        
          
    def fileDialog(self):
        self.directory = filedialog.askdirectory()  
        self.browseEntry.delete(first=0, last=1000)
        self.browseEntry.insert(tk.END, self.directory)        
  
    def accept(self):          
        # Get output directory        
        self.dir = self.browseEntry.get()
        if not self.dir:
            self.summaryDisplay.delete("1.0", tk.END)
            self.summaryDisplay.insert(tk.END, "Warning! You must select an output directory.\n\n")
            return
        # Get output filename
        self.name = self.fileEntry.get()
        if not self.name:
            self.summaryDisplay.delete("1.0", tk.END)
            self.summaryDisplay.insert(tk.END, "Warning! You must select an output NetCDF file name.\n\n")
            return
        if self.name[-3:] != ".nc":
            self.name = self.name + ".nc"        
        # Get complete filename
        self.name = self.dir + "/" + self.name
        # Get z-levels
        self.userz = self.zEntry.get().split(",")
        self.z = []
        for item in self.userz:
            try:
                self.z.append(-abs(float(item)))
            except ValueError:
                pass 
        self.z.sort(reverse=True)
        # Get user-selected boundaries
        self.northBoundary = int(self.northEntry.get())
        self.southBoundary = int(self.southEntry.get())
        self.eastBoundary  = int(self.eastEntry.get())
        self.westBoundary  = int(self.westEntry.get())
        # Parse variables
        self.userdict = {}
        for element in self.keylist:
            if element in self.userkeys:
                settings = list(self.vardict[element])
                settings[2] = self.descriptions[self.userkeys.index(element)]
                settings[7] = self.scales[self.userkeys.index(element)]
                settings[8] = [self.mincolran[self.userkeys.index(element)], \
                               self.maxcolran[self.userkeys.index(element)]]
                settings[9] = self.colors[self.userkeys.index(element)]
                self.userdict[element] = tuple(settings)
        # Get printing format 
        fmt = [None, "eps", "pdf", "png"][self.fmt.get()]
        
        self.set_widgets("disabled")
        
        if ( self.mode.get() < 2 ):
            status = self.layers(self.name, self.z, self.t0, self.t1, \
                            steply=self.steplySwitch.get(), daily=self.dailySwitch.get(), \
                            weekly=self.weeklySwitch.get(), monthly=self.monthlySwitch.get(), \
                            imonth=self.imonth.get(), iyear=self.iyear.get(), \
                            emonth=self.emonth.get(), eyear=self.eyear.get(), \
                            nb=self.northBoundary, sb=self.southBoundary, \
                            eb=self.eastBoundary,  wb=self.westBoundary, \
                            makeplot=self.makeplot.get(), fmt=fmt, userkey=self.userdict)
        else:  
            status = self.timeseries(self.name, self.z, self.t0, self.t1, \
                            steply=self.steplySwitch.get(), daily=self.dailySwitch.get(), \
                            weekly=self.weeklySwitch.get(), monthly=self.monthlySwitch.get(), \
                            imonth=self.imonth.get(), iyear=self.iyear.get(), \
                            emonth=self.emonth.get(), eyear=self.eyear.get(), \
                            nb=self.northBoundary, sb=self.southBoundary, \
                            eb=self.eastBoundary,  wb=self.westBoundary, \
                            makeplot=self.makeplot.get(), fmt=fmt, userkey=self.userdict)
                    
        if status:
            self.destroy()
        else:
            self.set_widgets("normal")                
         
    def set_widgets(self, mode):
        for w in self.browseFrame.winfo_children():
            w.configure(state=mode)
        for w in self.modeFrame.winfo_children():
            w.configure(state=mode)
        for w in self.quickFrame.winfo_children():
            w.configure(state=mode)                                
        for w in self.areaFrame.winfo_children():
            w.configure(state=mode)
        for w in self.manualFrame.winfo_children():
            try:
                w.configure(state=mode)
            except:
                pass
        for w in self.varFrame.winfo_children():
            try:
                w.configure(state=mode)
            except:
                pass
        for w in self.plotFrame.winfo_children():
            w.configure(state=mode)
        self.zEntry["state"] = mode
        self.ok["state"] = mode
        
    def varStruct(self, choices=()):
        from itertools import islice
        var = {
        "mean free surface":                  (("y", "x", "T"),       3, "mean sea surface height",                   "meter",         lambda x: x.mean(axis=2), ("zeta",),                 ("mean_free_surface",),             "linear", [-2, 2],      "jet", "real"), \
        "minimum free surface":               (("y", "x", "T"),       3, "minimum sea surface height",                "meter",         lambda x: x.min(axis=2),  ("zeta",),                 ("minimum_free_surface",),          "linear", [-2, 2],      "jet", "real"), \
        "maximum free surface":               (("y", "x", "T"),       3, "maximum sea surface height",                "meter",         lambda x: x.max(axis=2),  ("zeta",),                 ("maximum_free_surface",),          "linear", [-2, 2],      "jet", "real"), \
                    
        "average temperature":                (("y", "x", "z", "T"),  4, "mean potential temperature",                "Celsius",       lambda x: x.mean(axis=3), ("temp", "zeta"),          ("average_temperature",),           "linear", [5, 20],      "jet", "real"), \
    	"minimum temperature":                (("y", "x", "z", "T"),  4, "minimum potential temperature",             "Celsius",       lambda x: x.min(axis=3),  ("temp", "zeta"),          ("minimum_temperature",),           "linear", [5, 20],      "jet", "real"), \
    	"maximum temperature":                (("y", "x", "z", "T"),  4, "maximum potential temperature",             "Celsius",       lambda x: x.max(axis=3),  ("temp", "zeta"),          ("maximum_temperature",),           "linear", [5, 20],      "jet", "real"), \
        
    	"average surface temperature":        (("y", "x", "T"),       3, "mean surface temperature",                  "Celsius",       lambda x: x.mean(axis=2), ("temp",),                 ("average_surface_temperature",),   "linear", [10, 20],     "jet", "real"), \
    	"minimum surface temperature":        (("y", "x", "T"),       3, "minimum surface temperature",               "Celsius",       lambda x: x.min(axis=2),  ("temp",),                 ("minimum_surface_temperature",),   "linear", [10, 20],     "jet", "real"), \
    	"maximum surface temperature":        (("y", "x", "T"),       3, "maximum surface temperature",               "Celsius",       lambda x: x.max(axis=2),  ("temp",),                 ("maximum_surface_temperature",),   "linear", [10, 20],     "jet", "real"), \
        
    	"average bottom temperature":         (("y", "x", "T"),       3, "mean bottom potential temperature",         "Celsius",       lambda x: x.mean(axis=2), ("temp",),                 ("average_bottom_temperature",),    "linear", [2, 12],      "jet", "real"), \
    	"minimum bottom temperature":         (("y", "x", "T"),       3, "minimum bottom potential temperature",      "Celsius",       lambda x: x.min(axis=2),  ("temp",),                 ("minimum_bottom_temperature",),    "linear", [2, 12],      "jet", "real"), \
    	"maximum bottom temperature":         (("y", "x", "T"),       3, "maximum bottom potential temperature",      "Celsius",       lambda x: x.max(axis=2),  ("temp",),                 ("maximum_bottom_temperature",),    "linear", [2, 12],      "jet", "real"), \
        
        
    	"average salinity":                   (("y", "x", "z", "T"),  4, "mean salinity",                             "",              lambda x: x.mean(axis=3), ("salt", "zeta"),          ("average_salinity",),              "linear", [34.5, 35.5], "jet", "real"), \
    	"minimum salinity":                   (("y", "x", "z", "T"),  4, "minimum salinity",                          "",              lambda x: x.min(axis=3),  ("salt", "zeta"),          ("minimum_salinity",),              "linear", [34.5, 35.5], "jet", "real"), \
    	"maximum salinity":                   (("y", "x", "z", "T"),  4, "maximum salinity",                          "",              lambda x: x.max(axis=3),  ("salt", "zeta"),          ("maximum_salinity",),              "linear", [34.5, 35.5], "jet", "real"), \
        
    	"average surface salinity":           (("y", "x", "T"),       3, "mean surface salinity",                     "",              lambda x: x.mean(axis=2), ("salt",),                 ("average_surface_salinity",),      "linear", [34.5, 35.5], "jet", "real"), \
    	"minimum surface salinity":           (("y", "x", "T"),       3, "minimum surface salinity",                  "",              lambda x: x.min(axis=2),  ("salt",),                 ("minimum_surface_salinity",),      "linear", [34.5, 35.5], "jet", "real"), \
    	"maximum surface salinity":           (("y", "x", "T"),       3, "maximum surface salinity",                  "",              lambda x: x.max(axis=2),  ("salt",),                 ("maximum_surface_salinity",),      "linear", [34.5, 35.5], "jet", "real"), \
        
    	"average bottom salinity":            (("y", "x", "T"),       3, "mean bottom salinity",                      "",              lambda x: x.mean(axis=2), ("salt",),                 ("average_bottom_salinity",),       "linear", [34.5, 35.5], "jet", "real"), \
    	"minimum bottom salinity":            (("y", "x", "T"),       3, "minimum bottom salinity",                   "",              lambda x: x.min(axis=2),  ("salt",),                 ("minimum_bottom_salinity",),       "linear", [34.5, 35.5], "jet", "real"), \
    	"maximum bottom salinity":            (("y", "x", "T"),       3, "maximum bottom salinity",                   "",              lambda x: x.max(axis=2),  ("salt",),                 ("maximum_bottom_salinity",),       "linear", [34.5, 35.5], "jet", "real"), \
        
        
    	"average density":                    (("y", "x", "z", "T"),  4, "EOS-80 mean density",                       "kg m-3",        lambda x: x.mean(axis=3), ("temp", "salt", "zeta"),  ("average_density",),               "linear", [24, 29],     "jet", "real"), \
    	"minimum density":                    (("y", "x", "z", "T"),  4, "EOS-80 minimum density",                    "kg m-3",        lambda x: x.min(axis=3),  ("temp", "salt", "zeta"),  ("minimum_density",),               "linear", [24, 29],     "jet", "real"), \
    	"maximum density":                    (("y", "x", "z", "T"),  4, "EOS-80 maximum density",                    "kg m-3",        lambda x: x.max(axis=3),  ("temp", "salt", "zeta"),  ("maximum_density",),               "linear", [24, 29],     "jet", "real"), \
        
    	"average surface density":            (("y", "x", "T"),       3, "EOS-80 mean surface density",               "kg m-3",        lambda x: x.mean(axis=2), ("temp", "salt", "zeta"),  ("average_surface_density",),       "linear", [24, 29],     "jet", "real"), \
    	"minimum surface density":            (("y", "x", "T"),       3, "EOS-80 minimum surface density",            "kg m-3",        lambda x: x.min(axis=2),  ("temp", "salt", "zeta"),  ("minimum_surface_density",),       "linear", [24, 29],     "jet", "real"), \
    	"maximum surface density":            (("y", "x", "T"),       3, "EOS-80 maximum surface density",            "kg m-3",        lambda x: x.max(axis=2),  ("temp", "salt", "zeta"),  ("maximum_surface_density",),       "linear", [24, 29],     "jet", "real"), \
        
    	"average bottom density":             (("y", "x", "T"),       3, "EOS-80 mean bottom density",                "kg m-3",        lambda x: x.mean(axis=2), ("temp", "salt", "zeta"),  ("average_bottom_density",),        "linear", [24, 29],     "jet", "real"), \
    	"minimum bottom density":             (("y", "x", "T"),       3, "EOS-80 minimum bottom density",             "kg m-3",        lambda x: x.min(axis=2),  ("temp", "salt", "zeta"),  ("minimum_bottom_density",),        "linear", [24, 29],     "jet", "real"), \
    	"maximum bottom density":             (("y", "x", "T"),       3, "EOS-80 maximum bottom density",             "kg m-3",        lambda x: x.max(axis=2),  ("temp", "salt", "zeta"),  ("maximum_bottom_density",),        "linear", [24, 29],     "jet", "real"), \
        
        
        "average u":                          (("y", "x", "z", "T"),  4, "u-component of circulation (average)",      "cm s-1",        lambda x: x.mean(axis=3), ("u", "v", "zeta"),        ("average_u",),                     "linear", [-50, 50],    "bwr", "real"), \
	    "minimum u":                          (("y", "x", "z", "T"),  4, "u-component of circulation (minimum)",      "cm s-1",        lambda x: x.min(axis=3),  ("u", "v", "zeta"),        ("minimum_u",),                     "linear", [-50, 50],    "bwr", "real"), \
        "maximum u":                          (("y", "x", "z", "T"),  4, "u-component of circulation (maximum)",      "cm s-1",        lambda x: x.max(axis=3),  ("u", "v", "zeta"),        ("maximum_u",),                     "linear", [-50, 50],    "bwr", "real"), \
        
        "average surface u":                  (("y", "x", "T"),       3, "surface u-component (average)",             "cm s-1",        lambda x: x.mean(axis=2), ("u", "v"),                ("average_surface_u",),             "linear", [-50, 50],    "bwr", "real"), \
        "minimum surface u":                  (("y", "x", "T"),       3, "surface u-component (minimum)",             "cm s-1",        lambda x: x.min(axis=2),  ("u", "v"),                ("minimum_surface_u",),             "linear", [-50, 50],    "bwr", "real"), \
        "maximum surface u":                  (("y", "x", "T"),       3, "surface u-component (maximum)",             "cm s-1",        lambda x: x.max(axis=2),  ("u", "v"),                ("maximum_surface_u",),             "linear", [-50, 50],    "bwr", "real"), \
 
        "average bottom u":                   (("y", "x", "T"),       3, "bottom u-component (average)",              "cm s-1",        lambda x: x.mean(axis=2), ("u", "v"),                ("average_bottom_u",),              "linear", [-50, 50],    "bwr", "real"), \
        "minimum bottom u":                   (("y", "x", "T"),       3, "bottom u-component (minimum)",              "cm s-1",        lambda x: x.min(axis=2),  ("u", "v"),                ("minimum_bottom_u",),              "linear", [-50, 50],    "bwr", "real"), \
        "maximum bottom u":                   (("y", "x", "T"),       3, "bottom u-component (maximum)",              "cm s-1",        lambda x: x.max(axis=2),  ("u", "v"),                ("maximum_bottom_u",),              "linear", [-50, 50],    "bwr", "real"), \
	
    
        "average v":                          (("y", "x", "z", "T"),  4, "v-component of circulation (average)",      "cm s-1",        lambda x: x.mean(axis=3), ("u", "v", "zeta"),        ("average_v",),                     "linear", [-50, 50],    "bwr", "real"), \
	    "minimum v":                          (("y", "x", "z", "T"),  4, "v-component of circulation (minimum)",      "cm s-1",        lambda x: x.min(axis=3),  ("u", "v", "zeta"),        ("minimum_v",),                     "linear", [-50, 50],    "bwr", "real"), \
        "maximum v":                          (("y", "x", "z", "T"),  4, "v-component of circulation (maximum)",      "cm s-1",        lambda x: x.max(axis=3),  ("u", "v", "zeta"),        ("maximum_v",),                     "linear", [-50, 50],    "bwr", "real"), \
        
        "average surface v":                  (("y", "x", "T"),       3, "surface v-component (average)",             "cm s-1",        lambda x: x.mean(axis=2), ("u", "v"),                ("average_surface_v",),             "linear", [-50, 50],    "bwr", "real"), \
        "minimum surface v":                  (("y", "x", "T"),       3, "surface v-component (minimum)",             "cm s-1",        lambda x: x.min(axis=2),  ("u", "v"),                ("minimum_surface_v",),             "linear", [-50, 50],    "bwr", "real"), \
        "maximum surface v":                  (("y", "x", "T"),       3, "surface v-component (maximum)",             "cm s-1",        lambda x: x.max(axis=2),  ("u", "v"),                ("maximum_surface_v",),             "linear", [-50, 50],    "bwr", "real"), \
 
        "average bottom v":                   (("y", "x", "T"),       3, "bottom v-component (average)",              "cm s-1",        lambda x: x.mean(axis=2), ("u", "v"),                ("average_bottom_v",),              "linear", [-50, 50],    "bwr", "real"), \
        "minimum bottom v":                   (("y", "x", "T"),       3, "bottom v-component (minimum)",              "cm s-1",        lambda x: x.min(axis=2),  ("u", "v"),                ("minimum_bottom_v",),              "linear", [-50, 50],    "bwr", "real"), \
        "maximum bottom v":                   (("y", "x", "T"),       3, "bottom v-component (maximum)",              "cm s-1",        lambda x: x.max(axis=2),  ("u", "v"),                ("maximum_bottom_v",),              "linear", [-50, 50],    "bwr", "real"), \
                
        
        "average velocity":                      (("y", "x", "z", "T"),  4, "mean circulation velocity",                    "cm s-1",        lambda x: x.mean(axis=3), ("u", "v", "zeta"),  ("average_velocity_u", "average_velocity_v"),                 "linear", [0, 50],      "jet", "complex"), \
    	"maximum velocity":                      (("y", "x", "z", "T"),  4, "maximum circulation velocity",                 "cm s-1",        lambda x: x.max(axis=3),  ("u", "v", "zeta"),  ("maximum_magnitude", "maximum_direction"),                   "linear", [0, 50],      "jet", "complex"), \
            
        "average surface velocity":              (("y", "x", "T"),       3, "mean surface circulation velocity",            "cm s-1",        lambda x: x.mean(axis=2), ("u", "v"),          ("average_surface_velocity_u", "average_surface_velocity_v"), "linear", [0, 50],      "jet", "complex"), \
        "maximum surface velocity":              (("y", "x", "T"),       3, "maximum surface circulation velocity",         "cm s-1",        lambda x: x.max(axis=2),  ("u", "v"),          ("maximum_surface_magnitude", "maximum_surface_direction"),   "linear", [0, 50],      "jet", "complex"), \
    
        "average bottom velocity":               (("y", "x", "T"),       3, "mean bottom circulation velocity",             "cm s-1",        lambda x: x.mean(axis=2), ("u", "v"),          ("average_bottom_velocity_u", "average_bottom_velocity_v"),   "linear", [0, 50],      "jet", "complex"), \
        "maximum bottom velocity":               (("y", "x", "T"),       3, "maximum bottom circulation velocity",          "cm s-1",        lambda x: x.max(axis=2),  ("u", "v"),          ("maximum_bottom_magnitude", "maximum_bottom_direction"),     "linear", [0, 50],      "jet", "complex"), \
    
    
    	"average potential energy deficit":   (("y", "x", "T"),       3, "mean Potential Energy Deficit",             "kg m-1 s-2",    lambda x: x.mean(axis=2), ("temp", "salt", "zeta"),  ("average_potential_energy_deficit",),                        "log",    [1, 200],     "jet", "real"), \
    	"minimum potential energy deficit":   (("y", "x", "T"),       3, "minimum Potential Energy Deficit",          "kg m-1 s-2",    lambda x: x.min(axis=2),  ("temp", "salt", "zeta"),  ("minimum_potential_energy_deficit",),                        "log",    [1, 200],     "jet", "real"), \
    	"maximum potential energy deficit":   (("y", "x", "T"),       3, "maximum Potential Energy Deficit",          "kg m-1 s-2",    lambda x: x.max(axis=2),  ("temp", "salt", "zeta"),  ("maximum_potential_energy_deficit",),                        "log",    [1, 200],     "jet", "real"), \
    
    
    	"average mixed layer depth":          (("y", "x", "T"),       3, "mean Mixed Layer Depth",                    "meter",         lambda x: x.mean(axis=2), ("temp", "zeta"),          ("average_mixed_layer_depth",),                               "linear", [-100, 0],    "jet", "real"), \
    	"deepest mixed layer depth":          (("y", "x", "T"),       3, "deepest Mixed Layer Depth",                 "meter",         lambda x: x.min(axis=2),  ("temp", "zeta"),          ("deepest_mixed_layer_depth",),                               "linear", [-100, 0],    "jet", "real"), \
    	"shallowest mixed layer depth":       (("y", "x", "T"),       3, "shallowest Mixed Layer Depth",              "meter",         lambda x: x.max(axis=2),  ("temp", "zeta"),          ("shallowest_mixed_layer_depth",),                            "linear", [-100, 0],    "jet", "real"), \
    
    
    	"SST-based front index (average)":    (("y", "x", "T"),       3, "SST - Front Index (average)",               "",   lambda x: x.mean(axis=2), ("temp",),                            ("average_front_index",),                                     "linear", [0, 5],       "jet", "real"), \
    	"SST-based front index (minimum)":    (("y", "x", "T"),       3, "SST - Front Index (minimum)",               "",   lambda x: x.min(axis=2),  ("temp",),                            ("minimum_front_index",),                                     "linear", [0, 5],       "jet", "real"), \
    	"SST-based front index (maximum)":    (("y", "x", "T"),       3, "SST - Front Index (maximum)",               "",   lambda x: x.max(axis=2),  ("temp",),                            ("maximum_front_index",),                                     "linear", [0, 5],       "jet", "real"), \
        }
            
        # Subset from user's choices
        new = {}
        if choices:        
            for i in choices:
                val = var[next(islice(var, i, None))]
                key = list(var.keys())[list(var.values()).index(val)]
                new[key] = val
        if new: var = new
        return var      
        
    def layers(self, cdf, z, T0, T1, \
               steply=0, daily=0, weekly=0, monthly=0, imonth="Jan", iyear="1997", emonth="Dec", eyear="2015", 
               nb=90, sb=-90, eb=180, wb=-180, makeplot=0, fmt = "png", userkey=[]):
        
        """  """
                
        from datetime import timedelta        
        import numpy as np
        import os
        
        """ Return boolean """
        boolean = True
    
        """ Length of time dimension """
        T = len(T0)        
        
        """ Set defaults """
        # Default output variables
        var = self.varStruct(self.choices)
        
        """ Processing input variables """
        if userkey:
            var = userkey
        
        """ Extract fields """
        # Keyword
        key = list(var)
        # Number of dimensions
        ndim = [item[1] for item in var.values()]
        # Title
        long = [item[2] for item in var.values()]
        # Units
        units = [item[3] for item in var.values()]
        # Method
        func = [item[4] for item in var.values()]
        # ROMS associated variables
        roms = list(dict.fromkeys([item for group in [item[5] for item in var.values()] for item in group]))
        # Array
        array = [item[6] for item in var.values()];       
        # Color scale
        colorScale = [item[7] for item in var.values()]
        # Color axis
        colorAxis = [item[8] for item in var.values()]
        # Color map
        colorMap = [item[9] for item in var.values()]
        # Type
        tipo = [item[10] for item in var.values()]
        
        """ Depth """
        z = list(dict.fromkeys([-abs(z) for z in z]))        
        # Number of user-selected z-levels
        N = len(z)
        
        """ Show summary """
        self.summaryDisplay.delete("1.0", tk.END)
        self.summaryDisplay.insert(tk.END, "Summary of actions to be taken:\n")
        self.summaryDisplay.insert(tk.END, "\n OUTPUT: An output file " + cdf + " will be created\n")
        self.summaryDisplay.insert(tk.END, "\n MODE: Time averages (2-D maps) will be produced\n")
        if z:
            self.summaryDisplay.insert(tk.END, "\n Z-SLICING: z-slices of 3D variables will be obtained at constant depths [meters]:  " \
                                   + ', '.join([str(elem) for elem in z]) + "\n")
        
        TIME = []            
        """ Processing input dates (manual selection) """        
        for i in range(T):
            time = []
            # Start date 
            idate = T0[i]
            # End date 
            edate = T1[i]
            if idate > edate: 
                idate, edate = edate, idate
            if idate == edate: 
                time.extend([idate])                
                TIME.append(time)
                continue            
            # Build time vector
            while idate <= edate:
                time.append(idate)
                idate += timedelta(seconds=self.step_config)            
            TIME.append(time)
                
        """ Processing input dates (quick selection) """
        if ( imonth and iyear and emonth and eyear ):
            # Process starting month and year
            imonth = self.monthsList.index(imonth) + 1; iyear = int(iyear)
            # Process end month and year
            emonth = self.monthsList.index(emonth) + 1; eyear = int(eyear)            
            # Get starting serial month number
            imn = iyear * 12 + imonth
            # Get end serial month number
            emn = eyear * 12 + emonth
            if imn > emn:
                imn, emn = emn, imn       
            if steply:
                idate, edate = self.get_date(imn, emn)
                for fecha in self.dates:
                    if ( fecha >= idate and fecha <= edate):
                        TIME.append([fecha])
            if daily:
                idate, edate = self.get_date(imn, emn)
                time = []
                while idate <= ( edate + timedelta(seconds=self.step_config)):  
                    past = idate - timedelta(seconds=self.step_config)
                    if ( past.day != idate.day ):
                        if len(time) == 86400 / self.step_config:
                            TIME.append(time)
                        time = []
                    time.append(idate)
                    idate += timedelta(seconds=self.step_config)
            if weekly:     
                idate, edate = self.get_date(imn, emn)                    
                time = []
                while idate <= ( edate + timedelta(seconds=self.step_config)):         
                    past = idate - timedelta(seconds=self.step_config)
                    if ( idate.weekday() == 6 and past.weekday() != 6 ):                                                 
                        if len(time) == 7 * 86400 / self.step_config:                            
                            TIME.append(time)
                        time = []
                    time.append(idate)
                    idate += timedelta(seconds=self.step_config)            
            if monthly:
                while imn <= emn:
                    time = []
                    idate, edate = self.get_date(imn, imn)  
                    # Build time vector
                    while idate <= edate:
                        time.append(idate)
                        idate += timedelta(seconds=self.step_config)
                    TIME.append(time)
                    imn += 1
        
        if len(TIME) < 1: 
            self.summaryDisplay.delete("1.0", tk.END)
            self.summaryDisplay.insert(tk.END, "Warning! You must select at least one time period.\n\n")
            return 0
                    
        """ Display processing information """
        self.summaryDisplay.insert(tk.END, "\n TIME: The following time intervals will be considered:\n")
        for time in TIME:
            if ( time[0] == time[-1] ):
                self.summaryDisplay.insert(tk.END, "        {} (snapshot)\n".format(time[0].strftime("%d-%b-%Y %H:%M")))
            else:
                s_idate = time[ 0].strftime("%d-%b-%Y %H:%M")
                s_edate = time[-1].strftime("%d-%b-%Y %H:%M")
                self.summaryDisplay.insert(tk.END, "   from {} to {}\n".format(s_idate, s_edate))
    
        """ Fix unsorted boundaries, if needed """
        if sb > nb: 
            nb, sb = sb, nb
        if wb > eb: 
            eb, wb = wb, eb
        self.summaryDisplay.insert(tk.END, "\n BOUNDARIES: Processing will be limited to the following boundaries:\n\n" + \
                                   "   NORTH " + str(nb) + " j-grid" + \
                                   "   SOUTH " + str(sb) + " j-grid" + \
                                   "   EAST "  + str(abs(eb)) + " i-grid" + \
                                   "   WEST "  + str(abs(wb)) + " i-grid\n")        
        nb += 1; eb += 1
        self.summaryDisplay.insert(tk.END, "\n VARIABLES: The following variables will be included:\n\n")
        for i1, i4, i5, i6 in zip(key, colorScale, colorMap, colorAxis):
            string = "   '" + i1 + "'. Plot (if chosen) using " + i4 + " scale, '"  \
            + i5 + "' map and axis range from " + str(i6[0]) + " to " + str(i6[1]) + "\n\n" 
            self.summaryDisplay.insert(tk.END, string)
        
        if makeplot:
            self.summaryDisplay.insert(tk.END, " PLOTTING: yes, using " + fmt + \
                                       " format. Figures will be saved in an IMAGES directory\n\n")                    
        else:
            self.summaryDisplay.insert(tk.END, " PLOTTING: no\n\n")        
    
        self.summaryDisplay.insert(tk.END, "Continue? (Yes/No)")
        self.wait_variable(self.usersInput)
        if not self.usersInput.get():
            boolean = False
            return boolean
        
        # Get start times 
        T0 = [item[0]  for item in TIME]
        # Get end times
        T1 = [item[-1] for item in TIME]                           
        
        """ Get input file name matching initial date """
        self.rep = {"%Y": TIME[0][0].strftime("%Y"), \
               "%y": TIME[0][0].strftime("%y"), \
               "%m": TIME[0][0].strftime("%m"), \
               "%b": TIME[0][0].strftime("%b"), \
               "%d": TIME[0][0].strftime("%d"), \
               "%o": str(TIME[0][0].timetuple().tm_yday).zfill(3), \
               "%H": TIME[0][0].strftime("%H"), \
               "%M": TIME[0][0].strftime("%M"), \
               "%S": TIME[0][0].strftime("%S"), \
               "%rt": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.tstep_config, \
                                             TIME[0][0]), \
               "%rm": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.mstep_config, \
                                             TIME[0][0])}
        self.rep = dict((re.escape(kk), vv) for kk, vv in self.rep.items()) 
        pattern = re.compile("|".join(self.rep.keys()))
        f = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.tracers_convention)
        f = netCDF4.Dataset(f, "r")        
           
        """ Read variables """
        # Longitude
        x = f.variables["lon_rho"][sb : nb, wb : eb]   
        # Latitude
        y = f.variables["lat_rho"][sb : nb, wb : eb]
        # S-RHO
        s = f.variables["s_rho"][:]
        # Bathymetry
        H = f.variables["h"][sb : nb, wb : eb].filled(fill_value=np.nan)
        # Mask
        mask = f.variables["mask_rho"][sb : nb, wb : eb]
        # Angle
        gridangle = f.variables["angle"][sb : nb, wb : eb]
                
        """ Vertical parameters """
        # Vertical transformation equation
        Vtransform = f.variables["Vtransform"][:]
        # Vertical stretching function
        Vstretching = f.variables["Vstretching"][:]
        # Surface vertical stretching parameter
        theta_s = f.variables["theta_s"][:]
        # Bottom vertical stretching parameter
        theta_b = f.variables["theta_b"][:]
        # Critical depth
        hc = f.variables["hc"][:]  
        f.close()
        
        """ Get dimensions """
        Lp, Mp = x.shape        
        # Number of vertical layers
        Np = len(s)
        # Get quiver spacing
        SP = round(.5 * (Lp + Mp) / 40)
        
        """ Make NetCDF """
        self.usersInput.set(-99)
        if os.path.isfile(cdf):
            
            self.summaryDisplay.delete("1.0", tk.END)
            self.summaryDisplay.insert(tk.END, "Warning! " + cdf + \
                                       " file is going to be overwritten\n\nContinue? (Yes/No)")
            self.wait_variable(self.usersInput)
            if not self.usersInput.get():
                boolean = False
                return boolean        
            
            os.remove(cdf)
        print("STEP 0/" + str(1 + makeplot) + ": Building NetCDF...\n")
        self.summaryDisplay.delete("1.0", tk.END)
        self.summaryDisplay.insert(tk.END, "STEP 0/" + str(1 + makeplot) + ": Building NetCDF...\n\n")
        self.update(); self.summaryDisplay.yview_moveto(1)
        self.makecdf(cdf, (nb, sb, eb, wb), Lp, Mp, N, x, y, z, T0, T1, self.offset, var)
        o = netCDF4.Dataset(cdf, "a")
        
        """ MAIN LOOP """   
        print("\n\nSTEP 1/" + str(1 + makeplot) + ": Loop over time...\n")
        self.summaryDisplay.insert(tk.END, "\n\nSTEP 1/" + str(1 + makeplot) + ": Loop over time...\n\n")
        self.update(); self.summaryDisplay.yview_moveto(1)
        for item in TIME:
            print("   Processing from " + \
                           item[0].strftime("%d-%b-%Y %H:%M") + " to " + \
                           item[-1].strftime("%d-%b-%Y %H:%M") + "...")
            self.summaryDisplay.insert(tk.END, "   Processing from " + \
                           item[0].strftime("%d-%b-%Y %H:%M") + " to " + \
                           item[-1].strftime("%d-%b-%Y %H:%M") + "...\n")
            self.update(); self.summaryDisplay.yview_moveto(1)
            
            """ Initalize output arrays """
            for en in key:
                if en == "mean free surface":
                    mean_free_surface = []
                elif en == "minimum free surface":
                    minimum_free_surface = []                
                elif en == "maximum free surface":
                    maximum_free_surface = []                    
                elif en == "average temperature":                    
                    average_temperature = []                         
                elif en == "minimum temperature":
                    minimum_temperature = []
                elif en == "maximum temperature":
                    maximum_temperature = []                
                elif en == "average surface temperature":
                    average_surface_temperature = []                    
                elif en == "minimum surface temperature":
                    minimum_surface_temperature = []                    
                elif en == "maximum surface temperature":
                    maximum_surface_temperature = []                
                elif en == "average bottom temperature":
                    average_bottom_temperature = []
                elif en == "minimum bottom temperature":
                    minimum_bottom_temperature = []                
                elif en == "maximum bottom temperature":
                    maximum_bottom_temperature = []                
                elif en == "average salinity":
                    average_salinity = []                
                elif en == "minimum salinity":
                    minimum_salinity = []
                elif en == "maximum salinity":
                    maximum_salinity = []
                elif en == "average surface salinity":
                    average_surface_salinity = []
                elif en == "minimum surface salinity":
                    minimum_surface_salinity = []
                elif en == "maximum surface salinity":
                    maximum_surface_salinity = []                
                elif en == "average bottom salinity":
                    average_bottom_salinity = []                
                elif en == "minimum bottom salinity":
                    minimum_bottom_salinity = []                
                elif en == "maximum bottom salinity":
                    maximum_bottom_salinity = []                
                elif en == "average density":
                    average_density = []                
                elif en == "minimum density":
                    minimum_density = []                
                elif en == "maximum density":
                    maximum_density = []                
                elif en == "average surface density":
                    average_surface_density = []                
                elif en == "minimum surface density":
                    minimum_surface_density = []                
                elif en == "maximum surface density":
                    maximum_surface_density = []                
                elif en == "average bottom density":
                    average_bottom_density = []                
                elif en == "minimum bottom density":
                    minimum_bottom_density = []                
                elif en == "maximum bottom density":
                    maximum_bottom_density = []                
                elif en == "average u":
                    average_u = []                         
                elif en == "minimum u":
                    minimum_u = []                
                elif en == "maximum u":
                    maximum_u = []                         
                elif en == "average surface u":
                    average_surface_u = []                    
                elif en == "minimum surface u":
                    minimum_surface_u = []                
                elif en == "maximum surface u":
                    maximum_surface_u = []            
                elif en == "average bottom u":
                    average_bottom_u = []                    
                elif en == "minimum bottom u":
                    minimum_bottom_u = []            
                elif en == "maximum bottom u":
                    maximum_bottom_u = []                
                elif en == "average v":
                    average_v = []                         
                elif en == "minimum v":
                    minimum_v = []                
                elif en == "maximum v":
                    maximum_v = []            
                elif en == "average surface v":
                    average_surface_v = []                    
                elif en == "minimum surface v":
                    minimum_surface_v = []                
                elif en == "maximum surface v":
                    maximum_surface_v = []            
                elif en == "average bottom v":
                    average_bottom_v = []                    
                elif en == "minimum bottom v":
                    minimum_bottom_v = []            
                elif en == "maximum bottom v":
                    maximum_bottom_v = []  
                elif en == "average velocity":
                    average_velocity_u = []
                    average_velocity_v = []
                elif en == "maximum velocity":
                    maximum_magnitude = []            
                    maximum_direction = np.zeros((Lp, Mp, N))                
                elif en == "average surface velocity":
                    average_surface_velocity_u = []
                    average_surface_velocity_v = []
                elif en == "maximum surface velocity":
                    maximum_surface_magnitude = []
                    maximum_surface_direction = np.zeros((Lp, Mp))                    
                elif en == "average bottom velocity":
                    average_bottom_velocity_u = []
                    average_bottom_velocity_v = []
                elif en == "maximum bottom velocity":
                    maximum_bottom_magnitude = []  
                    maximum_bottom_direction = np.zeros((Lp, Mp))
                elif en == "average potential energy deficit":
                    average_potential_energy_deficit = []                
                elif en == "minimum potential energy deficit":
                    minimum_potential_energy_deficit = []                
                elif en == "maximum potential energy deficit":
                    maximum_potential_energy_deficit = []                
                elif en == "average mixed layer depth":
                    average_mixed_layer_depth = []
                elif en == "deepest mixed layer depth":
                    deepest_mixed_layer_depth = []                
                elif en == "shallowest mixed layer depth":
                    shallowest_mixed_layer_depth = []
                elif en == "SST-based front index (average)":
                    average_front_index = []                
                elif en == "SST-based front index (minimum)":
                    minimum_front_index = []                
                elif en == "SST-based front index (maximum)":
                    maximum_front_index = []                                    
                else:                              
                    raise ValueError("Unavailable variable")    
               
            COUNT = 0                
            for i in item:  
                print("\n" + i.strftime("%d-%b-%Y %H:%M"))
                self.summaryDisplay.insert(tk.END, "\n" + i.strftime("%d-%b-%Y %H:%M") + "\n")
                self.update(); self.summaryDisplay.yview_moveto(1)
                """ Get input file name matching date """
                self.rep = {"%Y": i.strftime("%Y"), \
                   "%y": i.strftime("%y"), \
                   "%m": i.strftime("%m"), \
                   "%b": i.strftime("%b"), \
                   "%d": i.strftime("%d"), \
                   "%o": str(i.timetuple().tm_yday).zfill(3), \
                   "%H": i.strftime("%H"), \
                   "%M": i.strftime("%M"), \
                   "%S": i.strftime("%S"), \
                   "%rt": self.get_running_number(self.todatetime(self.idate_config), \
                                                 self.tstep_config, i), \
                   "%rm": self.get_running_number(self.todatetime(self.idate_config), \
                                                 self.mstep_config, i)}
                self.rep = dict((re.escape(kk), vv) for kk, vv in self.rep.items()) 
                pattern = re.compile("|".join(self.rep.keys()))
                f = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.tracers_convention)                
                g = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.momentum_convention)                
                
                if ("temp" in roms or "salt" in roms):
                    try:
                        F = netCDF4.Dataset(f, "r")
                        print("   File " + f + " opened successfully for reading tracers.")
                        self.summaryDisplay.insert(tk.END, "   File " + f + " opened successfully for reading tracers.\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                    except FileNotFoundError:
                        print("   File " + f + " not found. Skipping...")
                        self.summaryDisplay.insert(tk.END, "   File " + f + " not found. Skipping...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)                    
                        continue
                if ("u" in roms or "v" in roms or "zeta" in roms):    
                    try:
                        G = netCDF4.Dataset(g, "r")
                        print("   File " + g + " opened successfully for reading momentum.")
                        self.summaryDisplay.insert(tk.END, "   File " + g + " opened successfully for reading momentum.\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)                    
                    except FileNotFoundError:
                        print("   File " + g + " not found. Skipping...")
                        self.summaryDisplay.insert(tk.END, "   File " + g + " not found. Skipping...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)                     
                        continue                
                
                COUNT += 1
                """ Dealing with time... """
                # Read time from ocean file
                if ("temp" in roms or "salt" in roms):
                    time_file = F.variables["ocean_time"][:].tolist()
                else:
                    time_file = G.variables["ocean_time"][:].tolist()
                # Convert current loop index (i) according to the offset in config file
                ioffset = (i - self.offset).total_seconds()
                # Find reading time index
                ti = time_file.index(ioffset)
                
                # Read variables
                for variable in roms:
                    if variable == "temp":
                        print("   Reading temperature...")
                        self.summaryDisplay.insert(tk.END, "   Reading temperature...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        temp = np.squeeze(F.variables["temp"][ti, :, sb : nb, wb : eb]).filled(fill_value=np.nan)
                    elif variable == "salt":
                        print("   Reading salinity...")
                        self.summaryDisplay.insert(tk.END, "   Reading salinity...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        salt = np.squeeze(F.variables["salt"][ti, :, sb : nb, wb : eb]).filled(fill_value=np.nan)
                    elif variable == "zeta":
                        print("   Reading free-surface...")
                        self.summaryDisplay.insert(tk.END, "   Reading free-surface...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        zeta = np.squeeze(G.variables["zeta"][ti, sb : nb, wb : eb]).filled(fill_value=np.nan)
                        # Compute depth at RHO points
                        print("   Computing depth at RHO-points...")
                        self.summaryDisplay.insert(tk.END, "   Computing depth at RHO-points...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        z_rho = self.zlevs(Vtransform, Vstretching, H, zeta, theta_s, theta_b, hc, Np, "r")    
                        # Compute depth at W points
                        print("   Computing depth at W-points...")
                        self.summaryDisplay.insert(tk.END, "   Computing depth at W-points...\n")   
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        z_w = self.zlevs(Vtransform, Vstretching, H, zeta, theta_s, theta_b, hc, Np, "w")
                    elif variable == "u":
                        print("   Reading u...")
                        self.summaryDisplay.insert(tk.END, "   Reading u...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        u = 100 * self.u2rho(np.squeeze(G.variables["u"][ti, :, :, :])).filled(fill_value=np.nan)
                        u = u[:, sb : nb, wb : eb]                    
                    elif variable == "v":
                        print("   Reading v...")
                        self.summaryDisplay.insert(tk.END, "   Reading v...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        v = 100 * self.v2rho(np.squeeze(G.variables["v"][ti, :, :, :])).filled(fill_value=np.nan)
                        v = v[:, sb : nb, wb : eb]
                
                if ("temp" in roms or "salt" in roms):
                    F.close()
                if ("u" in roms or "v" in roms or "zeta" in roms):  
                    G.close()
                 
                if ( "u" in roms or "v" in roms ):
                    # Take grid rotation into account
                    uang =  u * np.cos(gridangle) -  v * np.sin(gridangle)
                    vang =  u * np.sin(gridangle) +  v * np.cos(gridangle)
                    # Update u, v
                    u, v = uang, vang             
            
                for en in key:
                    if en == "mean free surface":
                        mean_free_surface = self.ProcessPeriod(mean_free_surface, \
                                                               zeta, COUNT, "avg")[0]                        
                    
                    elif en == "minimum free surface":
                        minimum_free_surface = self.ProcessPeriod(minimum_free_surface, \
                                                                  zeta, COUNT, "min")[0]
                    
                    elif en == "maximum free surface":
                        maximum_free_surface = self.ProcessPeriod(maximum_free_surface, \
                                                                  zeta, COUNT, "max")[0]
                        
                    elif en == "average temperature":                    
                        if "temp_z" not in locals():
                            print("   Temperature z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Temperature z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            temp_z = self.zslice(z_rho, temp, z, H, zeta) 
                        average_temperature = self.ProcessPeriod(average_temperature, \
                                                                 temp_z, COUNT, "avg")[0]                        
                             
                    elif en == "minimum temperature":
                        if "temp_z" not in locals():
                            print("   Temperature z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Temperature z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            temp_z = self.zslice(z_rho, temp, z, H, zeta)                    
                        minimum_temperature = self.ProcessPeriod(minimum_temperature, \
                                                                 temp_z, COUNT, "min")[0]
                    
                    elif en == "maximum temperature":
                        if "temp_z" not in locals():
                            print("   Temperature z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Temperature z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            temp_z = self.zslice(z_rho, temp, z, H, zeta)                    
                        maximum_temperature = self.ProcessPeriod(maximum_temperature, \
                                                                 temp_z, COUNT, "min")[0]
                    
                    elif en == "average surface temperature":
                        average_surface_temperature = self.ProcessPeriod(average_surface_temperature, \
                                                                         temp[-1, :, :], COUNT, "avg")[0]                       
                        
                    elif en == "minimum surface temperature":
                        minimum_surface_temperature = self.ProcessPeriod(minimum_surface_temperature, \
                                                                         temp[-1, :, :], COUNT, "min")[0]                        
                        
                    elif en == "maximum surface temperature":
                        maximum_surface_temperature = self.ProcessPeriod(maximum_surface_temperature, \
                                                                         temp[-1, :, :], COUNT, "max")[0]                        
                    
                    elif en == "average bottom temperature":
                        average_bottom_temperature = self.ProcessPeriod(average_bottom_temperature, \
                                                                        temp[0, :, :], COUNT, "avg")[0]
                    
                    elif en == "minimum bottom temperature":
                        minimum_bottom_temperature = self.ProcessPeriod(minimum_bottom_temperature, \
                                                                        temp[0, :, :], COUNT, "min")[0]
                    
                    elif en == "maximum bottom temperature":
                        maximum_bottom_temperature = self.ProcessPeriod(maximum_bottom_temperature, \
                                                                        temp[0, :, :], COUNT, "max")[0]
                    
                    elif en == "average salinity":                    
                        if "salt_z" not in locals():
                            print("   Salinity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Salinity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            salt_z = self.zslice(z_rho, salt, z, H, zeta) 
                        average_salinity = self.ProcessPeriod(average_salinity, \
                                                                 salt_z, COUNT, "avg")[0]                        
                             
                    elif en == "minimum salinity":
                        if "salt_z" not in locals():
                            print("   Salinity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Salinity z-slicing...\n")[0]
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            salt_z = self.zslice(z_rho, salt, z, H, zeta)                    
                        minimum_salinity = self.ProcessPeriod(minimum_salinity, \
                                                                 salt_z, COUNT, "min")[0]
                    
                    elif en == "maximum salinity":
                        if "salt_z" not in locals():
                            print("   Salinity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Salinity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            salt_z = self.zslice(z_rho, salt, z, H, zeta)                    
                        maximum_salinity = self.ProcessPeriod(maximum_salinity, \
                                                                 salt_z, COUNT, "min")[0]
                    
                    elif en == "average surface salinity":
                        average_surface_salinity = self.ProcessPeriod(average_surface_salinity, \
                                                                         salt[-1, :, :], COUNT, "avg")[0]                        
                        
                    elif en == "minimum surface salinity":
                        minimum_surface_salinity = self.ProcessPeriod(minimum_surface_salinity, \
                                                                         salt[-1, :, :], COUNT, "min")[0]                       
                        
                    elif en == "maximum surface salinity":
                        maximum_surface_salinity = self.ProcessPeriod(maximum_surface_salinity, \
                                                                         salt[-1, :, :], COUNT, "max")[0]                      
                    
                    elif en == "average bottom salinity":
                        average_bottom_salinity = self.ProcessPeriod(average_bottom_salinity, \
                                                                        salt[0, :, :], COUNT, "avg")[0]
                    
                    elif en == "minimum bottom salinity":
                        minimum_bottom_salinity = self.ProcessPeriod(minimum_bottom_salinity, \
                                                                        salt[0, :, :], COUNT, "min")[0]
                    
                    elif en == "maximum bottom salinity":
                        maximum_bottom_salinity = self.ProcessPeriod(maximum_bottom_salinity, \
                                                                        salt[0, :, :], COUNT, "max")[0]                    
                    
                    elif en == "average density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        if "dens_z" not in locals():
                            print("   Density z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Density z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens_z = self.zslice(z_rho, dens, z, H, zeta)
                        average_density = self.ProcessPeriod(average_density, \
                                                                 dens_z, COUNT, "avg")[0]  
                    
                    elif en == "minimum density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        if "dens_z" not in locals():
                            print("   Density z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Density z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens_z = self.zslice(z_rho, dens, z, H, zeta)
                        minimum_density = self.ProcessPeriod(minimum_density, \
                                                                 dens_z, COUNT, "min")[0]
                    
                    elif en == "maximum density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        if "dens_z" not in locals():
                            print("   Density z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   Density z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens_z = self.zslice(z_rho, dens, z, H, zeta)
                        maximum_density = self.ProcessPeriod(maximum_density, \
                                                                 dens_z, COUNT, "max")[0]                            
                                        
                    elif en == "average surface density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        average_surface_density = self.ProcessPeriod(average_surface_density, \
                                                                         dens[-1, :, :], COUNT, "avg")[0] 
                    
                    elif en == "minimum surface density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        minimum_surface_density = self.ProcessPeriod(minimum_surface_density, \
                                                                         dens[-1, :, :], COUNT, "min")[0] 
                    
                    elif en == "maximum surface density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        maximum_surface_density = self.ProcessPeriod(maximum_surface_density, \
                                                                         dens[-1, :, :], COUNT, "max")[0] 
                    
                    elif en == "average bottom density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        average_bottom_density = self.ProcessPeriod(average_bottom_density, \
                                                                         dens[0, :, :], COUNT, "avg")[0]                         
                  
                    
                    elif en == "minimum bottom density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        minimum_bottom_density = self.ProcessPeriod(minimum_bottom_density, \
                                                                         dens[0, :, :], COUNT, "min")[0]                                                
                    
                    elif en == "maximum bottom density":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        maximum_bottom_density = self.ProcessPeriod(maximum_bottom_density, \
                                                                         dens[0, :, :], COUNT, "max")[0]                            
                    
                    elif en == "average u":
                        if "u_z" not in locals():
                            print("   u-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            u_z = self.zslice(z_rho, u, z, H, zeta)                    
                        average_u = self.ProcessPeriod(average_u, \
                                                                 u_z, COUNT, "avg")[0]                          
                             
                    elif en == "minimum u":
                        if "u_z" not in locals():
                            print("   u-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            u_z = self.zslice(z_rho, u, z, H, zeta)      
                        minimum_u = self.ProcessPeriod(minimum_u, \
                                                                 u_z, COUNT, "min")[0]                                                    
                    
                    elif en == "maximum u":
                        if "u_z" not in locals():
                            print("   u-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            u_z = self.zslice(z_rho, u, z, H, zeta)    
                        maximum_u = self.ProcessPeriod(maximum_u, \
                                                                 u_z, COUNT, "max")[0]
                             
                    elif en == "average surface u":
                        average_surface_u = self.ProcessPeriod(average_surface_u, \
                                                                         u[-1, :, :], COUNT, "avg")[0]                         
                        
                    elif en == "minimum surface u":
                        minimum_surface_u = self.ProcessPeriod(minimum_surface_u, \
                                                                         u[-1, :, :], COUNT, "min")[0]                         
                                            
                    elif en == "maximum surface u":
                        maximum_surface_u = self.ProcessPeriod(maximum_surface_u, \
                                                                         u[-1, :, :], COUNT, "max")[0]                          
                
                    elif en == "average bottom u":                        
                        average_bottom_u = self.ProcessPeriod(average_bottom_u, \
                                                                         u[0, :, :], COUNT, "avg")[0]
                        
                    elif en == "minimum bottom u":
                        minimum_bottom_u = self.ProcessPeriod(minimum_bottom_u, \
                                                                         u[0, :, :], COUNT, "min")[0]                                           
                
                    elif en == "maximum bottom u":
                        maximum_bottom_u = self.ProcessPeriod(maximum_bottom_u, \
                                                                         u[0, :, :], COUNT, "max")[0]                                                                   
                    
                    elif en == "average v":
                        if "v_z" not in locals():
                            print("   v-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            v_z = self.zslice(z_rho, v, z, H, zeta)                    
                        average_v = self.ProcessPeriod(average_v, \
                                                                 v_z, COUNT, "avg")[0]                          
                             
                    elif en == "minimum v":
                        if "v_z" not in locals():
                            print("   v-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            v_z = self.zslice(z_rho, v, z, H, zeta)      
                        minimum_v = self.ProcessPeriod(minimum_v, \
                                                                 v_z, COUNT, "min")[0]                                                    
                    
                    elif en == "maximum v":
                        if "v_z" not in locals():
                            print("   v-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            v_z = self.zslice(z_rho, v, z, H, zeta)    
                        maximum_v = self.ProcessPeriod(maximum_v, \
                                                                 v_z, COUNT, "max")[0]
                             
                    elif en == "average surface v":
                        average_surface_v = self.ProcessPeriod(average_surface_v, \
                                                                         v[-1, :, :], COUNT, "avg")[0]                         
                        
                    elif en == "minimum surface v":
                        minimum_surface_v = self.ProcessPeriod(minimum_surface_v, \
                                                                         v[-1, :, :], COUNT, "min")[0]                          
                                            
                    elif en == "maximum surface v":
                        maximum_surface_v = self.ProcessPeriod(maximum_surface_v, \
                                                                         v[-1, :, :], COUNT, "max")[0]                          
                
                    elif en == "average bottom v":                        
                        average_bottom_v = self.ProcessPeriod(average_bottom_v, \
                                                                         v[0, :, :], COUNT, "avg")[0]
                        
                    elif en == "minimum bottom v":
                        minimum_bottom_v = self.ProcessPeriod(minimum_bottom_v, \
                                                                         v[0, :, :], COUNT, "min")[0]                                           
                
                    elif en == "maximum bottom v":
                        maximum_bottom_v = self.ProcessPeriod(maximum_bottom_v, \
                                                                         v[0, :, :], COUNT, "max")[0]                    
                        
                    elif en == "average velocity":
                        if "u_z" not in locals():
                            print("   u-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            u_z = self.zslice(z_rho, u, z, H, zeta)
                        average_velocity_u = self.ProcessPeriod(average_velocity_u, \
                                                                 u_z, COUNT, "avg")[0] 
                        if "v_z" not in locals():
                            print("   v-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            v_z = self.zslice(z_rho, v, z, H, zeta)   
                        average_velocity_v = self.ProcessPeriod(average_velocity_v, \
                                                                 v_z, COUNT, "avg")[0] 
                    
                    elif en == "maximum velocity":
                        if "u_z" not in locals():
                            print("   u-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            u_z = self.zslice(z_rho, u, z, H, zeta)
                        if "v_z" not in locals():
                            print("   v-velocity z-slicing...")
                            self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            v_z = self.zslice(z_rho, v, z, H, zeta)                                           
                        maximum_magnitude, index = self.ProcessPeriod(maximum_magnitude, \
                                                  (u_z**2 + v_z**2)**.5, COUNT, "max")
                        maximum_direction[index] = np.degrees(np.arctan2(v_z, u_z))[index] 
                        
                
                    elif en == "average surface velocity":
                        average_surface_velocity_u = self.ProcessPeriod(average_surface_velocity_u, \
                                                                         u[-1, :, :], COUNT, "avg")[0]
                        average_surface_velocity_v = self.ProcessPeriod(average_surface_velocity_v, \
                                                                         v[-1, :, :], COUNT, "avg")[0]        
        
                    elif en == "maximum surface velocity":
                        maximum_surface_magnitude, index = self.ProcessPeriod(maximum_surface_magnitude, \
                                                  (u[-1, :, :]**2 + v[-1, :, :]**2)**.5, COUNT, "max")
                        maximum_surface_direction[index] = np.degrees(np.arctan2(v[-1, :, :], u[-1, :, :]))[index] 
    
                    elif en == "average bottom velocity":
                        average_bottom_velocity_u = self.ProcessPeriod(average_bottom_velocity_u, \
                                                                         u[0, :, :], COUNT, "avg")[0]
                        average_bottom_velocity_v = self.ProcessPeriod(average_bottom_velocity_v, \
                                                                         v[0, :, :], COUNT, "avg")[0]  
            
                    elif en == "maximum bottom velocity":
                        maximum_bottom_magnitude, index = self.ProcessPeriod(maximum_bottom_magnitude, \
                                                  (u[0, :, :]**2 + v[0, :, :]**2)**.5, COUNT, "max")
                        maximum_bottom_direction[index] = np.degrees(np.arctan2(v[0, :, :], u[0, :, :]))[index]                         

                    
                    elif en == "average potential energy deficit":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        if "PED" not in locals():
                            print("   Computing PED...")
                            self.summaryDisplay.insert(tk.END, "   Computing PED...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            PED = self.potential_energy_deficit(dens, z_rho, z_w, zeta, 200)
                        average_potential_energy_deficit = self.ProcessPeriod(average_potential_energy_deficit, \
                                                                         PED, COUNT, "avg")[0]                         
                    
                    elif en == "minimum potential energy deficit":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        if "PED" not in locals():
                            print("   Computing PED...")
                            self.summaryDisplay.insert(tk.END, "   Computing PED...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            PED = self.potential_energy_deficit(dens, z_rho, z_w, zeta, 200)
                        minimum_potential_energy_deficit = self.ProcessPeriod(minimum_potential_energy_deficit, \
                                                                         PED, COUNT, "min")[0] 
                    
                    elif en == "maximum potential energy deficit":
                        if "dens" not in locals():
                            print("   Computing density...")
                            self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            dens = self.rho_eos(temp, salt, -z_rho) - 1000
                        if "PED" not in locals():
                            print("   Computing PED...")
                            self.summaryDisplay.insert(tk.END, "   Computing PED...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            PED = self.potential_energy_deficit(dens, z_rho, z_w, zeta, 200)
                        maximum_potential_energy_deficit = self.ProcessPeriod(maximum_potential_energy_deficit, \
                                                                         PED, COUNT, "max")[0] 
                    
                    elif en == "average mixed layer depth":
                        if "MLD" not in locals():
                            print("   Computing MLD...")
                            self.summaryDisplay.insert(tk.END, "   Computing MLD...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            MLD = self.mixed_layer_depth(z_rho, temp, H, mask)
                        average_mixed_layer_depth = self.ProcessPeriod(average_mixed_layer_depth, \
                                                                         MLD, COUNT, "avg")[0]                
                    
                    elif en == "deepest mixed layer depth":
                        if "MLD" not in locals():
                            print("   Computing MLD...")
                            self.summaryDisplay.insert(tk.END, "   Computing MLD...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            MLD = self.mixed_layer_depth(z_rho, temp, H, mask)
                        deepest_mixed_layer_depth = self.ProcessPeriod(deepest_mixed_layer_depth, \
                                                                         MLD, COUNT, "min")[0]       
                    
                    elif en == "shallowest mixed layer depth":
                        if "MLD" not in locals():
                            print("   Computing MLD...")
                            self.summaryDisplay.insert(tk.END, "   Computing MLD...\n")
                            self.update(); self.summaryDisplay.yview_moveto(1)
                            MLD = self.mixed_layer_depth(z_rho, temp, H, mask)
                        shallowest_mixed_layer_depth = self.ProcessPeriod(shallowest_mixed_layer_depth, \
                                                                         MLD, COUNT, "max")[0]    
                    
                    elif en == "SST-based front index (average)":
                        average_front_index = self.ProcessPeriod(average_front_index, \
                                                                 self.sobel(temp[-1, :, :], mask), COUNT, "avg")
                                            
                    elif en == "SST-based front index (minimum)":
                        minimum_front_index = self.ProcessPeriod(minimum_front_index, \
                                                                 self.sobel(temp[-1, :, :], mask), COUNT, "min")
                    
                    elif en == "SST-based front index (maximum)":
                        maximum_front_index = self.ProcessPeriod(maximum_front_index, \
                                                                 self.sobel(temp[-1, :, :], mask), COUNT, "max")
                        
                    else:                              
                        raise ValueError("Unavailable variable") 
                        
                # Remove temporal varibles from previous iteration
                if "temp_z" in locals(): del temp_z
                if "salt_z" in locals(): del salt_z
                if "dens"   in locals(): del dens
                if "dens_z" in locals(): del dens_z
                if "u_z"    in locals(): del u_z
                if "v_z"    in locals(): del v_z
                if "PED"    in locals(): del PED
                if "MLD"    in locals(): del MLD
                        
            """ Write output arrays """
            print("\n   Writing from " + \
                           item[0].strftime("%d-%b-%Y %H:%M") + " to " + \
                           item[-1].strftime("%d-%b-%Y %H:%M") + "...\n\n")
            self.summaryDisplay.insert(tk.END, "\n   Writing from " + \
                           item[0].strftime("%d-%b-%Y %H:%M") + " to " + \
                           item[-1].strftime("%d-%b-%Y %H:%M") + "...\n\n")
            self.update(); self.summaryDisplay.yview_moveto(1)            
            for en, nd, ar, ty in zip(key, ndim, array, tipo):
                if ty == "real":
                    if nd == 3:
                        o.variables[en][:, :, TIME.index(item)] = eval(ar[0])
                    elif nd ==4:
                        o.variables[en][:, :, :, TIME.index(item)] = eval(ar[0])
                else:
                    if nd == 3:
                        if en[0:3] == "ave":
                            U = eval(ar[0]); V = eval(ar[1])
                            # Write magnitude
                            o.variables[en][:, :, TIME.index(item), 0] = (U**2 + V**2)**.5
                            # Write direction
                            o.variables[en][:, :, TIME.index(item), 1] = np.degrees(np.arctan2(V, U))
                        elif en[0:3] == "max":
                            # Write magnitude
                            o.variables[en][:, :, TIME.index(item), 0] = eval(ar[0])
                            # Write direction
                            o.variables[en][:, :, TIME.index(item), 1] = eval(ar[1])
                    elif nd == 4:
                        if en[0:3] == "ave":
                            U = eval(ar[0]); V = eval(ar[1])
                            # Write magnitude
                            o.variables[en][:, :, :, TIME.index(item), 0] = (U**2 + V**2)**.5
                            # Write direction
                            o.variables[en][:, :, :, TIME.index(item), 1] = np.degrees(np.arctan2(V, U))
                        elif en[0:3] == "max":
                            # Write magnitude
                            o.variables[en][:, :, :, TIME.index(item), 0] = eval(ar[0])
                            # Write direction
                            o.variables[en][:, :, :, TIME.index(item), 1] = eval(ar[1])                            
                        
        o.close()
        o = netCDF4.Dataset(cdf, "r")
              
        if makeplot:
            import matplotlib.pyplot as plt
            import matplotlib.ticker as mticker 
            import matplotlib.colors as colors
            R = ( ( y.max() - y.min() ) / ( x.max() - x.min() ) ) # axis ratio
            # Create new directory
            if not os.path.isdir(self.dir + "/IMAGES"):
                os.mkdir(self.dir + "/IMAGES")
            print("\nSTEP 2/" + str(1 + makeplot) + ": Plotting...\n")
            self.summaryDisplay.insert(tk.END, "\n\nSTEP 2/" + str(1 + makeplot) + ": Plotting...\n\n")
            self.update(); self.summaryDisplay.yview_moveto(1)
            for en, nd, lg, un, cs, ca, mc, ty in zip(key, ndim, long, units, colorScale, colorAxis, colorMap, tipo):
                print("   " + en)
                self.summaryDisplay.insert(tk.END, "   " + en + "\n")
                self.update(); self.summaryDisplay.yview_moveto(1)
                var_i = o.variables[en][:]
                    
                if nd == 3:
                    for period in TIME:
                        # New axes
                        fig = plt.figure(figsize=(9.45, 9.45*R))
                        ax = fig.add_subplot(111)
                        ax.set_aspect("equal") 
                        
                        # x-axis
                        xtickslocs = [round(100*i)/100 for i in np.linspace(floor(x.min()), ceil(x.max()), 6)]
                        ax.xaxis.set_major_locator(mticker.FixedLocator(xtickslocs))                        
                        ax.set_xlim(left=x.min(), right=x.max()) 
                        ax.set_xticklabels([( str(abs(i)) + "°W" ) for i in xtickslocs], fontsize=12)                        
                        ax.xaxis.labelpad = 15
                        
                        # y-axis
                        ytickslocs = [round(100*i)/100 for i in np.linspace(floor(y.min()), ceil(y.max()), 6)]
                        ax.yaxis.set_major_locator(mticker.FixedLocator(ytickslocs))                                              
                        ax.set_ylim(bottom=y.min(), top=y.max())
                        ax.set_yticklabels([( str(i) + "°N" ) for i in ytickslocs], fontsize=12)
                        
                        # Draw coastline                    
                        self.drawshore(ax, self.coast_x, self.coast_y)
                        
                        # Draw colormap
                        if ( ty == "complex" ):                            
                            # Get magnitude
                            magn = var_i[:, :, TIME.index(period), 0]
                            # Get direction
                            dire = var_i[:, :, TIME.index(period), 1]
                            # Get (u, v)
                            u = np.cos(np.deg2rad(dire))
                            v = np.sin(np.deg2rad(dire))  
                            # Draw quiver plot
                            ax.quiver(x[::SP,::SP], y[::SP,::SP], u[::SP,::SP], v[::SP,::SP], scale=40, zorder=1)                            
                        else:
                            # Get magnitude
                            magn = var_i[:, :, TIME.index(period)]
                        
                        if cs == "linear":
                            h = ax.pcolor(x, y, magn, cmap=mc, vmin=ca[0], vmax=ca[1], zorder=0)
                        else:
                            h = ax.pcolor(x, y, magn, cmap=mc, norm=colors.LogNorm(vmin=ca[0], vmax=ca[1]), zorder=0)    
                                
                        # Set colorbar
                        cb = fig.colorbar(h)
                        cb.set_label(un, fontsize=16)
                        # cb.set_ticks(colorTicks)
                        cb.ax.tick_params(labelsize=16)
                        
                        # Set title                         
                        ax.set_title(lg, fontsize=12, pad=15)
                        if (period[0] == period[-1]):
                            ax.set_xlabel(period[0].strftime("%d-%b-%Y %H:%M"), fontsize=12)
                        else:                                     
                            ax.set_xlabel("from  " + period[0].strftime("%d-%b-%Y %H:%M") + \
                                      "  to  " + period[-1].strftime("%d-%b-%Y %H:%M"))
                            
                        plt.tight_layout()
                        # Save figure
                        plt.savefig(self.dir + "./IMAGES/" + en + "_" + period[0].strftime("%Y%m%d%H") + \
                                    "_" +  period[-1].strftime("%Y%m%d%H") + "." + fmt, dpi=144)
                        
                        # Close figure
                        plt.close()
                     
                        
                elif nd == 4:
                    for period in TIME:
                        for level in z:
                            # New axes
                            fig = plt.figure(figsize=(9.45, 9.45*R))
                            ax = fig.add_subplot(111)
                            ax.set_aspect("equal") 
                        
                            # x-axis
                            xtickslocs = [round(100*i)/100 for i in np.linspace(floor(x.min()), ceil(x.max()), 6)]
                            ax.xaxis.set_major_locator(mticker.FixedLocator(xtickslocs))                            
                            ax.set_xlim(left=x.min(), right=x.max()) 
                            ax.set_xticklabels([( str(abs(i)) + "°W" ) for i in xtickslocs], fontsize=12)  
                            ax.xaxis.labelpad = 15
                        
                            # y-axis
                            ytickslocs = [round(100*i)/100 for i in np.linspace(floor(y.min()), ceil(y.max()), 6)]
                            ax.yaxis.set_major_locator(mticker.FixedLocator(ytickslocs))                            
                            ax.set_ylim(bottom=y.min(), top=y.max())
                            ax.set_yticklabels([( str(i) + "°W" ) for i in ytickslocs], fontsize=12)                    
                            
                            # Draw coastline                    
                            self.drawshore(ax, self.coast_x, self.coast_y)
                            
                            if ( ty == "complex" ):                                
                                # Get magnitude
                                magn = var_i[:, :, z.index(level), TIME.index(period), 0]
                                # Get direction
                                dire = var_i[:, :, z.index(level), TIME.index(period), 1]
                                # Get (u, v)
                                u = np.cos(np.deg2rad(dire))
                                v = np.sin(np.deg2rad(dire))  
                                # Draw quiver plot
                                ax.quiver(x[::SP,::SP], y[::SP,::SP], u[::SP,::SP], v[::SP,::SP], scale=40, zorder=1)
                                
                            else:
                                # Get magnitude
                                magn = var_i[:, :, z.index(level), TIME.index(period)]
                        
                            if cs == "linear":
                                h = ax.pcolor(x, y, magn, cmap=mc, vmin=ca[0], vmax=ca[1], zorder=0)
                            else:
                                h = ax.pcolor(x, y, magn, cmap=mc, norm=colors.LogNorm(vmin=ca[0], vmax=ca[1]), zorder=0)                                                               
                                    
                            # Set colorbar
                            cb = fig.colorbar(h)
                            cb.set_label(un, fontsize=16)
                            # cb.set_ticks(colorTicks)
                            cb.ax.tick_params(labelsize=16)
                            
                            # Set title 
                            ax.set_title(lg, fontsize=12, pad=15)
                            if (period[0] == period[-1]):
                                ax.set_xlabel(period[0].strftime("%d-%b-%Y %H:%M") + " Z = " + str(int(-level)) + " m", fontsize=12)
                            else:
                                ax.set_xlabel("from  " + period[0].strftime("%d-%b-%Y %H:%M") + \
                                          "  to  " + period[-1].strftime("%d-%b-%Y %H:%M") + \
                                          " Z = " + str(int(-level)) + " m", fontsize=12)
                                
                            plt.tight_layout()
                            # Save figure
                            plt.savefig(self.dir + "./IMAGES/" + en + "_" + period[0].strftime("%Y%m%d%H") + \
                                        "_" +  period[-1].strftime("%Y%m%d%H") + "_Z = " + str(int(-level)) + "." + fmt, dpi=144)
                            
                            # Close figure
                            plt.close()                
        o.close()
        return boolean
    
    def timeseries(self, cdf, z, T0, T1, \
               steply=0, daily=0, weekly=0, monthly=0, imonth="Jan", iyear="1997", emonth="Dec", eyear="2015", 
               nb=90, sb=-90, eb=180, wb=-180, makeplot=0, fmt = "png", userkey=[]):        
        """  """
                
        from datetime import timedelta        
        import numpy as np
        import os
    
        """ Return boolean """
        boolean = True
    
        """ Length of time dimension """
        T = len(T0)        
        
        """ Set defaults """
        # Default output variables
        var = self.varStruct(self.choices)
        
        """ Processing input variables """
        if userkey:
            var = userkey
        
        """ Extract fields """
        # Keyword
        key = list(var)
        # Number of dimensions
        ndim = [item[1] for item in var.values()]
        # Title
        long = [item[2] for item in var.values()]
        # Units
        units = [item[3] for item in var.values()]
        # Method
        func = [item[4] for item in var.values()]
        # ROMS associated variables
        roms = list(dict.fromkeys([item for group in [item[5] for item in var.values()] for item in group]))
        # Array
        array = [item[6] for item in var.values()]; uarray = list(dict.fromkeys(array))  
        # Color scale
        colorScale = [item[7] for item in var.values()]
        # Color axis
        colorAxis = [item[8] for item in var.values()]
        # Color map
        colorMap = [item[9] for item in var.values()]
        # Type
        tipo = [item[10] for item in var.values()]
        
        """ Depth """
        z = list(dict.fromkeys([-abs(z) for z in z]))        
        # Number of user-selected z-levels
        N = len(z)
        
        """ Show summary """
        self.summaryDisplay.delete("1.0", tk.END)
        self.summaryDisplay.insert(tk.END, "Summary of actions to be taken:\n")
        self.summaryDisplay.insert(tk.END, "\n OUTPUT: An output file " + cdf + " will be created\n")
        self.summaryDisplay.insert(tk.END, "\n MODE: Time series will be produced\n")
        if z:
            self.summaryDisplay.insert(tk.END, "\n Z-SLICING: z-slices of 3D variables will be obtained at constant depths [meters]:  " \
                                   + ', '.join([str(elem) for elem in z]) + "\n")
        
        TIME = []            
        """ Processing input dates (manual selection) """        
        for i in range(T):
            time = []
            # Start date 
            idate = T0[i]
            # End date 
            edate = T1[i]
            if idate > edate: 
                idate, edate = edate, idate
            if idate == edate: 
                time.extend([idate])                
                TIME.append(time)
                continue            
            # Build time vector
            while idate <= edate:
                time.append(idate)
                idate += timedelta(seconds=self.step_config)            
            TIME.append(time)
                
        """ Processing input dates (quick selection) """
        if ( imonth and iyear and emonth and eyear ):
            # Process starting month and year
            imonth = self.monthsList.index(imonth) + 1; iyear = int(iyear)
            # Process end month and year
            emonth = self.monthsList.index(emonth) + 1; eyear = int(eyear)            
            # Get starting serial month number
            imn = iyear * 12 + imonth
            # Get end serial month number
            emn = eyear * 12 + emonth
            if imn > emn:
                imn, emn = emn, imn  
            if steply:
                idate, edate = self.get_date(imn, emn)
                for fecha in self.dates:
                    if ( fecha >= idate and fecha <= edate):
                        TIME.append([fecha])
            if daily:
                idate, edate = self.get_date(imn, emn)
                time = []
                while idate <= ( edate + timedelta(seconds=self.step_config)):  
                    past = idate - timedelta(seconds=self.step_config)
                    if ( past.day != idate.day ):
                        if len(time) == 86400 / self.step_config:
                            TIME.append(time)
                        time = []
                    time.append(idate)
                    idate += timedelta(seconds=self.step_config)                
            if weekly:     
                idate, edate = self.get_date(imn, emn)                    
                time = []
                while idate <= ( edate + timedelta(seconds=self.step_config)):         
                    past = idate - timedelta(seconds=self.step_config)
                    if ( idate.weekday() == 6 and past.weekday() != 6 ):                                                 
                        if len(time) == 7 * 86400 / self.step_config:                            
                            TIME.append(time)
                        time = []
                    time.append(idate)
                    idate += timedelta(seconds=self.step_config)            
            if monthly:
                while imn <= emn:
                    time = []
                    idate, edate = self.get_date(imn, imn)  
                    # Build time vector
                    while idate <= edate:
                        time.append(idate)
                        idate += timedelta(seconds=self.step_config)
                    TIME.append(time)
                    imn += 1
        
        if len(TIME) < 1: 
            self.summaryDisplay.delete("1.0", tk.END)
            self.summaryDisplay.insert(tk.END, "Warning! You must select at least one time period.\n\n")
            return 0
                    
        """ Display processing information """
        self.summaryDisplay.insert(tk.END, "\n TIME: The following time intervals will be considered:\n")
        for time in TIME:
            if ( time[0] == time[-1] ):
                self.summaryDisplay.insert(tk.END, "        {} (snapshot)\n".format(time[0].strftime("%d-%b-%Y %H:%M")))
            else:
                s_idate = time[ 0].strftime("%d-%b-%Y %H:%M")
                s_edate = time[-1].strftime("%d-%b-%Y %H:%M")
                self.summaryDisplay.insert(tk.END, "   from {} to {}\n".format(s_idate, s_edate))
    
        """ Fix unsorted boundaries, if needed """
        if sb > nb: 
            nb, sb = sb, nb
        if wb > eb: 
            eb, wb = wb, eb
        self.summaryDisplay.insert(tk.END, "\n BOUNDARIES: Processing will be limited to the following boundaries:\n\n" + \
                                   "   NORTH " + str(nb) + " j-grid" + \
                                   "   SOUTH " + str(sb) + " j-grid" + \
                                   "   EAST "  + str(abs(eb)) + " i-grid" + \
                                   "   WEST "  + str(abs(wb)) + " igrid\n")        
        nb += 1; eb += 1
        self.summaryDisplay.insert(tk.END, "\n VARIABLES: The following variables will be included:\n\n")
        for i1, i4, i5, i6 in zip(key, colorScale, colorMap, colorAxis):
            string = "   '" + i1 + "'. Plot (if chosen) using " + i4 + " scale, '"  \
            + i5 + "' map and axis range from " + str(i6[0]) + " to " + str(i6[1]) + "\n\n" 
            self.summaryDisplay.insert(tk.END, string)
        
        if makeplot:
            self.summaryDisplay.insert(tk.END, " PLOTTING: yes, using " + fmt + \
                                       " format. Figures will be saved in an IMAGES directory\n\n")                    
        else:
            self.summaryDisplay.insert(tk.END, " PLOTTING: no\n\n")        
    
        self.summaryDisplay.insert(tk.END, "Continue? (Yes/No)")
        self.wait_variable(self.usersInput)
        if not self.usersInput.get():
            boolean = False
            return boolean
        
        # Get start times 
        T0 = [item[0]  for item in TIME]
        # Get end times
        T1 = [item[-1] for item in TIME]                           
        # Convert into single list 
        time = [item for sublist in TIME for item in sublist]
        # Remove duplicates
        time = list(dict.fromkeys(time))
        # Sort
        time.sort()
        # Get number of iterations over time
        T = len(time)           
        
        """ Get input file name matching initial date """
        self.rep = {"%Y": time[0].strftime("%Y"), \
               "%y": time[0].strftime("%y"), \
               "%m": time[0].strftime("%m"), \
               "%b": time[0].strftime("%b"), \
               "%d": time[0].strftime("%d"), \
               "%o": str(time[0].timetuple().tm_yday).zfill(3), \
               "%H": time[0].strftime("%H"), \
               "%M": time[0].strftime("%M"), \
               "%S": time[0].strftime("%S"), \
               "%rt": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.tstep_config, \
                                             time[0]), \
               "%rm": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.mstep_config, \
                                             time[0])}
        self.rep = dict((re.escape(kk), vv) for kk, vv in self.rep.items()) 
        pattern = re.compile("|".join(self.rep.keys()))
        f = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.tracers_convention)
        f = netCDF4.Dataset(f, "r")  
           
        """ Read variables """
        # Longitude
        x = f.variables["lon_rho"][sb : nb, wb : eb]   
        # Latitude
        y = f.variables["lat_rho"][sb : nb, wb : eb]
        # S-RHO
        s = f.variables["s_rho"][:]
        # Bathymetry
        H = f.variables["h"][sb : nb, wb : eb].filled(fill_value=np.nan)
        # Mask
        mask = f.variables["mask_rho"][sb : nb, wb : eb]
        # Angle
        gridangle = f.variables["angle"][sb : nb, wb : eb]        
                
        """ Vertical parameters """
        # Vertical transformation equation
        Vtransform = f.variables["Vtransform"][:]
        # Vertical stretching function
        Vstretching = f.variables["Vstretching"][:]
        # Surface vertical stretching parameter
        theta_s = f.variables["theta_s"][:]
        # Bottom vertical stretching parameter
        theta_b = f.variables["theta_b"][:]
        # Critical depth
        hc = f.variables["hc"][:]
        
        """ Get dimensions """
        Lp, Mp = x.shape        
        # Number of vertical layers
        Np = len(s)
        
        """ Make NetCDF """
        self.usersInput.set(-99)
        if os.path.isfile(cdf):
            
            self.summaryDisplay.delete("1.0", tk.END)
            self.summaryDisplay.insert(tk.END, "Warning! " + cdf + \
                                       " file is going to be overwritten\n\nContinue? (Yes/No)")
            self.wait_variable(self.usersInput)
            if not self.usersInput.get():
                boolean = False
                return boolean        
            
            os.remove(cdf)
        print("STEP 0/" + str(1 + makeplot) + ": Building NetCDF...\n")
        self.summaryDisplay.delete("1.0", tk.END)
        self.summaryDisplay.insert(tk.END, "STEP 0/" + str(1 + makeplot) + ": Building NetCDF...\n\n")
        self.update(); self.summaryDisplay.yview_moveto(1)
        self.cdfseries(cdf, (nb, sb, eb, wb), N, T, z, time, self.offset, var)        
        o = netCDF4.Dataset(cdf, "a")        
        
        f.close()  
          
        """ MAIN LOOP """     
        print("\nSTEP 1/" + str(1 + makeplot) + ": Loop over time...\n")
        self.summaryDisplay.insert(tk.END, "\n\nSTEP 1/" + str(2 + makeplot) + ": Loop over time...\n\n")
        self.update()
        for item in time:
            print("\n" + i.strftime("%d-%b-%Y %H:%M"))            
            self.summaryDisplay.insert(tk.END, "\n" + item.strftime("%d-%b-%Y %H:%M") + "\n")
            self.update(); self.summaryDisplay.yview_moveto(1)
            
            """ Get input file name matching date """
            self.rep = {"%Y": item.strftime("%Y"), \
               "%y": item.strftime("%y"), \
               "%m": item.strftime("%m"), \
               "%b": item.strftime("%b"), \
               "%d": item.strftime("%d"), \
               "%o": str(item.timetuple().tm_yday).zfill(3), \
               "%H": item.strftime("%H"), \
               "%M": item.strftime("%M"), \
               "%S": item.strftime("%S"), \
               "%rt": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.tstep_config, \
                                             item), \
               "%rm": self.get_running_number(self.todatetime(self.idate_config), \
                                             self.mstep_config, \
                                             item)}
            self.rep = dict((re.escape(kk), vv) for kk, vv in self.rep.items()) 
            pattern = re.compile("|".join(self.rep.keys()))
            f = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.tracers_convention)
            g = pattern.sub(lambda m: self.rep[re.escape(m.group(0))], self.momentum_convention)
            
            if ("temp" in roms or "salt" in roms):
                try:
                    F = netCDF4.Dataset(f, "r")
                    print("   File " + f + " opened successfully for reading tracers.")
                    self.summaryDisplay.insert(tk.END, "   File " + f + " opened successfully for reading tracers.\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                except FileNotFoundError:
                    print("   File " + f + " not found. Skipping...")
                    self.summaryDisplay.insert(tk.END, "   File " + f + " not found. Skipping...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)                    
                    continue
            if ("u" in roms or "v" in roms or "zeta" in roms):
                try:
                    G = netCDF4.Dataset(g, "r")
                    print("   File " + g + " opened successfully for reading momentum.")
                    self.summaryDisplay.insert(tk.END, "   File " + g + " opened successfully for reading momentum.\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)                    
                except FileNotFoundError:
                    print("   File " + g + " not found. Skipping...")
                    self.summaryDisplay.insert(tk.END, "   File " + g + " not found. Skipping...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)                     
                    continue   
                
            """ Dealing with time... """
            # Read time from ocean file
            if ("temp" in roms or "salt" in roms):
                time_file = F.variables["ocean_time"][:].tolist()
            else:
                time_file = G.variables["ocean_time"][:].tolist()
            # Convert current loop index (item) according to the offset in config file
            ioffset = (item - self.offset).total_seconds()
            # Find reading time index
            ti = time_file.index(ioffset)
            
            # Read variables
            for variable in roms:
                if variable == "temp":
                    print("   Reading temperature...")
                    self.summaryDisplay.insert(tk.END, "   Reading temperature...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    temp = np.squeeze(F.variables["temp"][ti, :, sb : nb, wb : eb]).filled(fill_value=np.nan)
                elif variable == "salt":
                    print("   Reading salinity...")
                    self.summaryDisplay.insert(tk.END, "   Reading salinity...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    salt = np.squeeze(F.variables["salt"][ti, :, sb : nb, wb : eb]).filled(fill_value=np.nan)
                elif variable == "zeta":
                    print("   Reading free-surface...")
                    self.summaryDisplay.insert(tk.END, "   Reading free-surface...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    zeta = np.squeeze(G.variables["zeta"][ti, sb : nb, wb : eb]).filled(fill_value=np.nan)
                    # Compute depth at RHO points
                    print("   Computing depth at RHO-points...")
                    self.summaryDisplay.insert(tk.END, "   Computing depth at RHO-points...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    z_rho = self.zlevs(Vtransform, Vstretching, H, zeta, theta_s, theta_b, hc, Np, "r")    
                    # Compute depth at W points
                    print("   Computing depth at W-points...")
                    self.summaryDisplay.insert(tk.END, "   Computing depth at W-points...\n")   
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    z_w = self.zlevs(Vtransform, Vstretching, H, zeta, theta_s, theta_b, hc, Np, "w")
                elif variable == "u":
                    print("   Reading u...")
                    self.summaryDisplay.insert(tk.END, "   Reading u...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    u = 100 * self.u2rho(np.squeeze(G.variables["u"][ti, :, :, :])).filled(fill_value=np.nan)
                    u = u[:, sb : nb, wb : eb]
                elif variable == "v":
                    print("   Reading v...")
                    self.summaryDisplay.insert(tk.END, "   Reading v...\n")
                    self.update(); self.summaryDisplay.yview_moveto(1)
                    v = 100 * self.v2rho(np.squeeze(G.variables["v"][ti, :, :, :])).filled(fill_value=np.nan)
                    v = v[:, sb : nb, wb : eb]
            
            if ( "u" in roms or "v" in roms ):
                # Take grid rotation into account
                uang =  u * np.cos(gridangle) -  v * np.sin(gridangle)
                vang =  u * np.sin(gridangle) +  v * np.cos(gridangle)
                # Update u, v
                u, v = uang, vang            
            
            # Compute variables
            for en in key:
                if en == "mean free surface":
                    o.variables[en][time.index(item)] = zeta.mean()
                
                elif en == "minimum free surface":
                    o.variables[en][time.index(item)] = zeta.min()
                
                elif en == "maximum free surface":
                    o.variables[en][time.index(item)] = zeta.max()
                    
                elif en == "average temperature":                    
                    if "temp_z" not in locals():
                        print("   Temperature z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Temperature z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        temp_z = self.zslice(z_rho, temp, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = temp_z[:, :, lvl].mean()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value
                         
                elif en == "minimum temperature":
                    if "temp_z" not in locals():
                        print("   Temperature z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Temperature z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        temp_z = self.zslice(z_rho, temp, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = temp_z[:, :, lvl].min()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value                    
                
                elif en == "maximum temperature":
                    if "temp_z" not in locals():
                        print("   Temperature z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Temperature z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        temp_z = self.zslice(z_rho, temp, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = temp_z[:, :, lvl].max()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value
                
                elif en == "average surface temperature":
                    o.variables[en][time.index(item)] = temp[-1, :, :].mean()
                    
                elif en == "minimum surface temperature":
                    o.variables[en][time.index(item)] = temp[-1, :, :].min()
                    
                elif en == "maximum surface temperature":
                    o.variables[en][time.index(item)] = temp[-1, :, :].max()
                
                elif en == "average bottom temperature":
                    o.variables[en][time.index(item)] = temp[0, :, :].mean()
                
                elif en == "minimum bottom temperature":
                    o.variables[en][time.index(item)] = temp[0, :, :].min()
                
                elif en == "maximum bottom temperature":
                    o.variables[en][time.index(item)] = temp[0, :, :].max()
                
                elif en == "average salinity":
                    if "salt_z" not in locals():
                        print("   Salinity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Salinity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        salt_z = self.zslice(z_rho, salt, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = salt_z[:, :, lvl].mean()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value
                
                elif en == "minimum salinity":
                    if "salt_z" not in locals():
                        print("   Salinity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Salinity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        salt_z = self.zslice(z_rho, salt, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = salt_z[:, :, lvl].min()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value
                
                elif en == "maximum salinity":
                    if "salt_z" not in locals():
                        print("   Salinity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Salinity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        salt_z = self.zslice(z_rho, salt, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = salt_z[:, :, lvl].max()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value
                
                elif en == "average surface salinity":
                    o.variables[en][time.index(item)] = salt[-1, :, :].mean()
                
                elif en == "minimum surface salinity":
                    o.variables[en][time.index(item)] = salt[-1, :, :].min()
                
                elif en == "maximum surface salinity":
                    o.variables[en][time.index(item)] = salt[-1, :, :].max()
                
                elif en == "average bottom salinity":
                    o.variables[en][time.index(item)] = salt[0, :, :].mean()
                
                elif en == "minimum bottom salinity":
                    o.variables[en][time.index(item)] = salt[0, :, :].min()
                
                elif en == "maximum bottom salinity":
                    o.variables[en][time.index(item)] = salt[0, :, :].max()
                
                elif en == "average density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    if "dens_z" not in locals():
                        print("   Density z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Density z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens_z = self.zslice(z_rho, dens, z, H, zeta)
                    for lvl in list(range(N)):
                        value = dens_z[:, :, lvl].mean()
                        if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                        o.variables[en][lvl, time.index(item)] = value                    
                
                elif en == "minimum density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    if "dens_z" not in locals():
                        print("   Density z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Density z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens_z = self.zslice(z_rho, dens, z, H, zeta)
                    for lvl in list(range(N)):
                        value = dens_z[:, :, lvl].min()
                        if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                        o.variables[en][lvl, time.index(item)] = value                    
                
                elif en == "maximum density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    if "dens_z" not in locals():
                        print("   Density z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   Density z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens_z = self.zslice(z_rho, dens, z, H, zeta)
                    for lvl in list(range(N)):
                        value = dens_z[:, :, lvl].max()
                        if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                        o.variables[en][lvl, time.index(item)] = value                    
                
                elif en == "average surface density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    o.variables[en][time.index(item)] = dens[-1, :, :].mean()
                
                elif en == "minimum surface density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    o.variables[en][time.index(item)] = dens[-1, :, :].min()
                
                elif en == "maximum surface density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    o.variables[en][time.index(item)] = dens[-1, :, :].max()                    
                
                elif en == "average bottom density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    o.variables[en][time.index(item)] = dens[0, :, :].mean()                   
                
                elif en == "minimum bottom density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    o.variables[en][time.index(item)] = dens[0, :, :].min()                    
                
                elif en == "maximum bottom density":
                    if "dens" not in locals():
                        print("   Computing density...")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    o.variables[en][time.index(item)] = dens[0 :, :].max()
                
                elif en == "average u":
                    if "u_z" not in locals():
                        print("   u-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        u_z = self.zslice(z_rho, u, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = u_z[:, :, lvl].mean()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value  
                         
                elif en == "minimum u":
                    if "u_z" not in locals():
                        print("   u-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        u_z = self.zslice(z_rho, u, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = u_z[:, :, lvl].min()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value                          
                
                elif en == "maximum u":
                    if "u_z" not in locals():
                        print("   u-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        u_z = self.zslice(z_rho, u, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = u_z[:, :, lvl].max()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value
                         
                elif en == "average surface u":
                    o.variables[en][time.index(item)] = u[-1, :, :].mean()
                    
                elif en == "minimum surface u":
                    o.variables[en][time.index(item)] = u[-1, :, :].min()                    
                
                elif en == "maximum surface u":
                    o.variables[en][time.index(item)] = u[-1, :, :].max()
            
                elif en == "average bottom u":
                    o.variables[en][time.index(item)] = u[0, :, :].mean()
                    
                elif en == "minimum bottom u":
                    o.variables[en][time.index(item)] = u[0, :, :].min()                    
            
                elif en == "maximum bottom u":
                    o.variables[en][time.index(item)] = u[0, :, :].max()
                
                elif en == "average v":
                    if "v_z" not in locals():
                        print("  v-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        v_z = self.zslice(z_rho, v, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = v_z[:, :, lvl].mean()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value   
                         
                elif en == "minimum v":
                    if "v_z" not in locals():
                        print("  v-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        v_z = self.zslice(z_rho, v, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = v_z[:, :, lvl].min()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value                           
                
                elif en == "maximum v":
                    if "v_z" not in locals():
                        print("  v-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        v_z = self.zslice(z_rho, v, z, H, zeta)                    
                    for lvl in list(range(N)):
                         value = v_z[:, :, lvl].max()
                         if type(value) == np.ma.core.MaskedConstant:
                             value = np.nan
                         o.variables[en][lvl, time.index(item)] = value                      
            
                elif en == "average surface v":
                    o.variables[en][time.index(item)] = v[-1, :, :].mean() 
                    
                elif en == "minimum surface v":
                    o.variables[en][time.index(item)] = v[-1, :, :].min()                     
                
                elif en == "maximum surface v":
                    o.variables[en][time.index(item)] = v[-1, :, :].max()     
            
                elif en == "average bottom v":
                    o.variables[en][time.index(item)] = v[0, :, :].mean()   
                    
                elif en == "minimum bottom v":
                    o.variables[en][time.index(item)] = v[0, :, :].min()                     
            
                elif en == "maximum bottom v":
                    o.variables[en][time.index(item)] = v[0, :, :].max()  
                    
                elif en == "average velocity":
                    if "u_z" not in locals():
                        print("  u-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        u_z = self.zslice(z_rho, u, z, H, zeta)
                    if "v_z" not in locals():
                        print("  v-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        v_z = self.zslice(z_rho, v, z, H, zeta)   
                    for lvl in list(range(N)):
                         u_mean = u_z[:, :, lvl].mean()
                         v_mean = v_z[:, :, lvl].mean()
                         w = (u_mean**2 + v_mean**2)**.5                         
                         if type(w) == np.ma.core.MaskedConstant:
                             w = np.nan
                         o.variables[en][lvl, time.index(item), 0] = w
                         ANG = np.degrees(np.arctan2(v_mean, u_mean))
                         if type(ANG) == np.ma.core.MaskedConstant:
                             ANG = np.nan
                         o.variables[en][lvl, time.index(item), 1] = ANG
                
                elif en == "maximum velocity":
                    if "u_z" not in locals():
                        print("  u-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   u-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        u_z = self.zslice(z_rho, u, z, H, zeta)
                    if "v_z" not in locals():
                        print("  v-velocity z-slicing...")
                        self.summaryDisplay.insert(tk.END, "   v-velocity z-slicing...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        v_z = self.zslice(z_rho, v, z, H, zeta)                      
                    for lvl in list(range(N)):
                        w = (u_z[:, :, lvl]**2 + v_z[:, :, lvl]**2)**.5
                        w_max = w.max()
                        if type(w_max) == np.ma.core.MaskedConstant:
                            w_max = np.nan
                        o.variables[en][lvl, time.index(item), 0] = w_max
                        ANG = np.degrees(np.arctan2(v_z[:, :, lvl], u_z[:, :, lvl]))
                        ANG = ANG[np.unravel_index(np.argmax(w), w.shape)]
                        if type(ANG) == np.ma.core.MaskedConstant:
                            ANG = np.nan
                        o.variables[en][lvl, time.index(item), 1] = ANG   
            
                elif en == "average surface velocity":
                    u_mean = u[-1, :, :].mean()
                    v_mean = v[-1, :, :].mean()
                    w = (u_mean**2 + v_mean**2)**.5     
                    if type(w) == np.ma.core.MaskedConstant:
                             w = np.nan
                    o.variables[en][time.index(item), 0] = w
                    ANG = np.degrees(np.arctan2(v_mean, u_mean))
                    if type(ANG) == np.ma.core.MaskedConstant:
                         ANG = np.nan
                    o.variables[en][time.index(item), 1] = ANG
    
                elif en == "maximum surface velocity":
                    w = (u[-1, :, :]**2 + v[-1, :, :]**2)**.5
                    w_max = w.max()
                    if type(w_max) == np.ma.core.MaskedConstant:
                        w_max = np.nan
                    o.variables[en][time.index(item), 0] = w_max
                    ANG = np.degrees(np.arctan2(v[-1, :, :], u[-1, :, :]))
                    ANG = ANG[np.unravel_index(np.argmax(w), w.shape)]
                    if type(ANG) == np.ma.core.MaskedConstant:
                        ANG = np.nan
                    o.variables[en][time.index(item), 1] = ANG 

                elif en == "average bottom velocity":
                    u_mean = u[0, :, :].mean()
                    v_mean = v[0, :, :].mean()
                    w = (u_mean**2 + v_mean**2)**.5     
                    if type(w) == np.ma.core.MaskedConstant:
                             w = np.nan
                    o.variables[en][time.index(item), 0] = w
                    ANG = np.degrees(np.arctan2(v_mean, u_mean))
                    if type(ANG) == np.ma.core.MaskedConstant:
                         ANG = np.nan
                    o.variables[en][time.index(item), 1] = ANG
        
                elif en == "maximum bottom velocity":
                    w = (u[0, :, :]**2 + v[0, :, :]**2)**.5
                    w_max = w.max()
                    if type(w_max) == np.ma.core.MaskedConstant:
                        w_max = np.nan
                    o.variables[en][time.index(item), 0] = w_max
                    ANG = np.degrees(np.arctan2(v[0, :, :], u[0, :, :]))
                    ANG = ANG[np.unravel_index(np.argmax(w), w.shape)]
                    if type(ANG) == np.ma.core.MaskedConstant:
                        ANG = np.nan
                    o.variables[en][time.index(item), 1] = ANG 
                
                elif en == "average potential energy deficit":
                    if "dens" not in locals():
                        print("   Computing density")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    if "PED" not in locals():
                        print("   Computing PED")
                        self.summaryDisplay.insert(tk.END, "   Computing PED...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        PED = self.potential_energy_deficit(dens, z_rho, z_w, zeta, 200)
                    o.variables[en][time.index(item)] = PED.mean() 
                
                elif en == "minimum potential energy deficit":
                    if "dens" not in locals():
                        print("   Computing density")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    if "PED" not in locals():
                        print("   Computing PED")
                        self.summaryDisplay.insert(tk.END, "   Computing PED...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        PED = self.potential_energy_deficit(dens, z_rho, z_w, zeta, 200)
                    o.variables[en][time.index(item)] = PED.min() 
                
                elif en == "maximum potential energy deficit":
                    if "dens" not in locals():
                        print("   Computing density")
                        self.summaryDisplay.insert(tk.END, "   Computing density...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        dens = self.rho_eos(temp, salt, -z_rho) - 1000
                    if "PED" not in locals():
                        print("   Computing PED")
                        self.summaryDisplay.insert(tk.END, "   Computing PED...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        PED = self.potential_energy_deficit(dens, z_rho, z_w, zeta, 200)
                    o.variables[en][time.index(item)] = PED.max() 
                
                elif en == "average mixed layer depth":
                    if "MLD" not in locals():
                        print("   Computing MLD")
                        self.summaryDisplay.insert(tk.END, "   Computing MLD...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        MLD = self.mixed_layer_depth(z_rho, temp, H, mask)
                    o.variables[en][time.index(item)] = MLD.mean()                 
                
                elif en == "deepest mixed layer depth":
                    if "MLD" not in locals():
                        print("   Computing MLD")
                        self.summaryDisplay.insert(tk.END, "   Computing MLD...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        MLD = self.mixed_layer_depth(z_rho, temp, H, mask)
                    o.variables[en][time.index(item)] = MLD.min()      
                
                elif en == "shallowest mixed layer depth":
                    if "MLD" not in locals():
                        print("   Computing MLD")
                        self.summaryDisplay.insert(tk.END, "   Computing MLD...\n")
                        self.update(); self.summaryDisplay.yview_moveto(1)
                        MLD = self.mixed_layer_depth(z_rho, temp, H, mask)
                    o.variables[en][time.index(item)] = MLD.max()    
                
                elif en == "SST-based front index (average)":
                    o.variables[en][time.index(item)] = self.sobel(temp[-1, :, :], mask).mean()
                
                elif en == "SST-based front index (minimum)":
                    o.variables[en][time.index(item)] = self.sobel(temp[-1, :, :], mask).min()
                
                elif en == "SST-based front index (maximum)":
                    o.variables[en][time.index(item)] = self.sobel(temp[-1, :, :], mask).max()                    
                
                else:                              
                    raise ValueError("Unavailable variable")  
                    
            # Remove temporal varibles from previous iteration
            if "temp_z" in locals(): del temp_z
            if "salt_z" in locals(): del salt_z
            if "dens"   in locals(): del dens
            if "dens_z" in locals(): del dens_z
            if "u_z"    in locals(): del u_z
            if "v_z"    in locals(): del v_z
            if "w"      in locals(): del w
            if "w_z"    in locals(): del w_z
            if "PED"    in locals(): del PED
            if "MLD"    in locals(): del MLD
                            
            f.close()    
            g.close()
                        
        o.close()
                
        if makeplot:
            
            o = netCDF4.Dataset(cdf, "r"); 
            
            tiempo = o.variables["time"][:]; tiempo = [self.offset + timedelta(seconds=item) for item in tiempo]
            
            from matplotlib.dates import  DateFormatter
            import matplotlib.pyplot as plt         
                        
            # Create new directory
            if not os.path.isdir(self.dir + "/IMAGES"):
                os.mkdir(self.dir + "/IMAGES")
                                   
            print("\nSTEP 2/" + str(1 + makeplot) + ": Plotting...\n") 
            self.summaryDisplay.insert(tk.END, "\n\nSTEP 2/" + str(1 + makeplot) + ": Plotting...\n\n")
            self.update(); self.summaryDisplay.yview_moveto(1)
            for en, nd, lg, un, cs, ca, mc, ty in zip(key, ndim, long, units, colorScale, colorAxis, colorMap, tipo):
                print("   " + en)
                self.summaryDisplay.insert(tk.END, "   " + en + "\n")
                self.update(); self.summaryDisplay.yview_moveto(1)
                var_i = o.variables[en][:].tolist()
               
                if nd == 3:
                    for period in TIME:
                        if ty == "complex":
                            data = [(i, var_i[tiempo.index(i)][0]) for i in tiempo if i in period]                            
                        else:
                            data = [(i, var_i[tiempo.index(i)]) for i in tiempo if i in period]
                        
                        # New axes
                        fig, ax = plt.subplots(figsize=(13.11, 8.10))                        
                        
                        # Plot time series
                        ax.plot_date(*zip(*data), "b.-"); ax.grid()
                        
                        if ty == "complex":                            
                            dire = [var_i[tiempo.index(i)][1] for i in tiempo if i in period]
                            u_plot = np.cos(np.deg2rad(dire)); v_plot = np.sin(np.deg2rad(dire))
                            ax.quiver(period, [0 for i in period], u_plot, v_plot)
                        
                        # x-axis
                        ax.xaxis.set_major_formatter( DateFormatter('%d-%b-%Y') ); plt.xticks(rotation=90)                                                
                        
                        # y-axis
                        plt.yscale(cs)
                        plt.ylim(ca[0], ca[1])
                        if cs == "linear":
                            ytickslocs = np.arange(ca[0], ca[1]+1e-3, (ca[1]-ca[0])/5)                        
                            plt.yticks(ticks=ytickslocs)
                        ax.set_ylabel(un, fontsize=16)
                                                                        
                        # Set title 
                        ax.set_title(lg, fontsize=12, pad=15)
                             
                        plt.tight_layout()
                        # Save figure
                        plt.savefig(self.dir + "./IMAGES/" + en + "_" + period[0].strftime("%Y%m%d%H") + \
                                    "_" +  period[-1].strftime("%Y%m%d%H") + "." + fmt, dpi=144)
                        
                        # Close figure
                        plt.close()
                        
                elif nd == 4:
                    for level in z:
                        if ty == "complex":
                            var_z = var_i[z.index(level)][:]
                        else:                                
                            var_z = var_i[z.index(level)]
                            
                        for period in TIME:
                            if ty == "complex":
                                data = [(i, var_z[tiempo.index(i)][0]) for i in tiempo if i in period]                            
                            else:
                                data = [(i, var_z[tiempo.index(i)]) for i in tiempo if i in period]
                            
                            # New axes
                            fig, ax = plt.subplots(figsize=(13.11, 8.10))                        
                            
                            # Plot time series
                            ax.plot_date(*zip(*data), "b.-"); ax.grid()
                            
                            if ty == "complex":                            
                                dire = [var_z[tiempo.index(i)][1] for i in tiempo if i in period]
                                u_plot = np.cos(np.deg2rad(dire)); v_plot = np.sin(np.deg2rad(dire))
                                ax.quiver(period, [0 for i in period], u_plot, v_plot)
                            
                            # x-axis
                            ax.xaxis.set_major_formatter( DateFormatter('%d-%b-%Y') ); plt.xticks(rotation=90) 
                                                        
                            # y-axis
                            plt.yscale(cs)
                            plt.ylim(ca[0], ca[1])
                            if cs == "linear":
                                ytickslocs = np.arange(ca[0], ca[1]+1e-3, (ca[1]-ca[0])/5)                        
                                plt.yticks(ticks=ytickslocs)
                            ax.set_ylabel(un, fontsize=16)
                                                                            
                            # Set title 
                            ax.set_title(lg + " Z = " + str(int(-level)) + " m", fontsize=12, pad=15)
                                
                            plt.tight_layout()
                            # Save figure
                            plt.savefig(self.dir + "./IMAGES/" + en + "_" + period[0].strftime("%Y%m%d%H") + \
                                        "_" +  period[-1].strftime("%Y%m%d%H") + "_Z = " + str(int(-level)) + "." + fmt, dpi=144)
                            
                            # Close figure
                            plt.close()                            
                            
                    
                                        
            o.close()
            
        return boolean
        
    def potential_energy_deficit(self, rho, z_rho, z_w, zeta, D):
        """ Potential Energy Deficit (Planque et al. 2006) """
        
        import numpy as np    
        # Get dimensions
        (N, M, L) = rho.shape
        
        z = np.tile((zeta - D), (N+1, 1, 1))
        # Get depth
        z_w[z_w < z] = z[z_w < z]
        # Compute layer thickness
        thick = z_w[1:, :, :] - z_w[:-1, :, :]
        # Get vertically-averaged water column density
        avg = (rho*thick).sum(axis=0)/thick.sum(axis=0)
        # Get Potential Energy Deficit
        PED = (9.80665*(avg-rho)*z_rho*thick).sum(axis=0)/thick.sum(axis=0)    
        return PED
    
    def mixed_layer_depth(self, z, temp, h, mask):
        """ Mixed Layer Depth: 
               
               First, the temperature difference between the water column and the 
               sea surface is calculated. The mixed layer depth is defined as the
               minimum depth where the temperature difference is 0.5 Celsius
               
            Arguments:
                
                Z      shape (N, M, L) depth (positive) at RHO points
                
                TEMP   shape (N, M, L) temperature
            
                H      shape (M, L) undisturbed water column depth
                
                MASK   shape (M, L) sea(1)/land(0) binary mask        
            """
        
        import numpy as np       
        # Subtract from sea surface temperature
        temp = temp[-1, :, :] - temp
        # Get dimensions 
        (N, M, L) = temp.shape
        # Output initialization
        mld = np.full((M, L), -99999.9)
        # Iterate over sigma layers, from bottom to top
        for i in range(N-1):
            # temperature below
            t0 = temp[i, :, :]
            # temperature above
            t1 = temp[i+1, :, :]
            # depth below
            z0 = z[i, :, :]
            # depth above
            z1 = z[i+1, :, :]
            # Check whether the MLD condition is met
            w = t0 >= .5
            # Vertical linear interpolation 
            mld[w] = z1[w] + (z0[w] - z1[w])*(.5 - t1[w])/(t0[w] - t1[w]) 
            
        # Where the MLD condition is not met...
        full_mixed = mld < -99999
        # ... set H as mixed layer depth
        mld[full_mixed] = -h[full_mixed]   
        # Mask
        mld[mask < .5] = 1e+6; mld[mld > 999999] = np.nan
        
        return mld
           
    def _dens0(self, S, T):
        """Density of seawater at zero pressure"""
    
        # --- Define constants ---
        a0 = 999.842594
        a1 =   6.793952e-2
        a2 =  -9.095290e-3
        a3 =   1.001685e-4
        a4 =  -1.120083e-6
        a5 =   6.536332e-9
    
        b0 =   8.24493e-1
        b1 =  -4.0899e-3
        b2 =   7.6438e-5
        b3 =  -8.2467e-7
        b4 =   5.3875e-9
    
        c0 =  -5.72466e-3
        c1 =   1.0227e-4
        c2 =  -1.6546e-6
    
        d0 =   4.8314e-4
    
        # --- Computations ---
        # Density of pure water
        SMOW = a0 + (a1 + (a2 + (a3 + (a4 + a5*T)*T)*T)*T)*T
    
        # More temperature polynomials
        RB = b0 + (b1 + (b2 + (b3 + b4*T)*T)*T)*T
        RC = c0 + (c1 + c2*T)*T
    
        return SMOW + RB*S + RC*(S**1.5) + d0*S*S
    
    # -----------------------------------------------------------------
    
    
    def _seck(self, S, T, P=0):
        """Secant bulk modulus"""
    
        # --- Pure water terms ---
    
        h0 =  3.239908
        h1 =  1.43713E-3
        h2 =  1.16092E-4
        h3 = -5.77905E-7
        AW = h0 + (h1 + (h2 + h3*T)*T)*T
    
        k0 =  8.50935E-5
        k1 = -6.12293E-6
        k2 =  5.2787E-8
        BW = k0 + (k1 + k2*T)*T
    
        e0 = 19652.21
        e1 = 148.4206
        e2 = -2.327105
        e3 =  1.360477E-2
        e4 = -5.155288E-5
        KW = e0 + (e1 + (e2 + (e3 + e4*T)*T)*T)*T
    
        # --- seawater, P = 0 ---
    
        SR = S**0.5
    
        i0 =  2.2838E-3
        i1 = -1.0981E-5
        i2 = -1.6078E-6
        j0 =  1.91075E-4
        A  = AW + (i0 + (i1 + i2*T)*T + j0*SR)*S
    
        f0 = 54.6746
        f1 = -0.603459
        f2 =  1.09987E-2
        f3 = -6.1670E-5
        g0 =  7.944E-2
        g1 =  1.6483E-2
        g2 = -5.3009E-4
        K0 = KW + (f0 + (f1 + (f2 + f3*T)*T)*T  \
                + (g0 + (g1 + g2*T)*T)*SR)*S
    
        # --- General expression ---
    
        m0 = -9.9348E-7
        m1 =  2.0816E-8
        m2 =  9.1697E-10
        B = BW + (m0 + (m1 + m2*T)*T)*S
    
        K = K0 + (A + B*P)*P
    
        return K
    
    # ----------------------------------------------
    
    def rho_eos(self, T, S, P=0):
        """Compute density of seawater from salinity, temperature, and pressure
        Usage: dens(S, T, [P])
        Input:
            S = Salinity,     [PSS-78]
            T = Temperature,  [Celsius]
            P = Pressure,     [dbar = 10**4 Pa]
        P is optional, with default value zero
        Output:
            Density,          [kg/m**3]
        Algorithm: UNESCO 1983
        """
    
        P = 0.1*P # Convert to bar
        return self._dens0(S,T)/(1 - P/self._seck(S,T,P))
        
    def sobel(self, X, mask):
        """ Sobel operator """
        import numpy 
        from scipy import ndimage
        
        dx = ndimage.sobel(X, 0)  
        dy = ndimage.sobel(X, 1)  
        
        grad = numpy.hypot(dx, dy)
        grad[mask < .5] = 1e+6; grad = numpy.ma.masked_array(grad, grad > 999999)
        
        grad[:,0] = numpy.nan; grad[:,-1] = numpy.nan; 
        grad[0,:] = numpy.nan; grad[-1,:] = numpy.nan;
        
        grad.mask[numpy.isnan(grad)] = True
        
        return grad
    
    def zslice(self, z, v, zq, H, zeta):
        import numpy as np
        # Get dimensions
        N, M, L = z.shape
        # Reshape
        z = np.reshape(z, (N, L * M))
        v = np.reshape(v, (N, L * M))
        H = -np.reshape(H, (L * M,))
        zeta = np.reshape(zeta, (L * M,))
        # Initialize output array
        vq = np.empty((len(zq), L * M)); vq[:] = np.nan
        
        for zi in zq:        
            # For each water column, find index of first depth greater than 'zi'
            gt = (z > zi).argmax(axis=0); 
            GT = gt + N * np.asarray(list(range(L * M)))
            if GT[0] < 1:
                GT[0] += 1
            # Find greater than 'zi'        
            zg = z[np.unravel_index(GT, z.shape, "F")]
            vg = v[np.unravel_index(GT, v.shape, "F")]
            # Find lower than 'zi'
            zl = z[np.unravel_index(GT-1, z.shape, "F")]
            vl = v[np.unravel_index(GT-1, v.shape, "F")]
            # Linearly interpolate to depth 'zi'
            vqi = vl + (vg - vl) * (zi - zl) / (zg - zl); vqi[gt == 0] = np.nan
            # If 'zi' is just above sea floor, get the deepest value
            vqi[np.logical_and(zi < z[0, :], zi >= H)] = v[0, np.logical_and(zi < z[0, :], zi >= H)]
            # If 'zi' is just below sea surface, get the shallowest value
            vqi[np.logical_and(zi > z[-1, :], zi <= zeta)] = v[-1, np.logical_and(zi > z[-1, :], zi <= zeta)]
            vq[zq.index(zi), :] = vqi
            
        # Reshape
        vq = np.transpose(np.reshape(vq, (len(zq), M, L)), (1, 2, 0)) 
        return vq
    
    def zlevs(self, Vtransform, Vstretching, H, zeta, theta_s, theta_b, hc, N, igrid):
        import numpy as np        
        # Get sigma
        if igrid == "w":
            sc = np.asarray([(k - N)/N      for k in range(0, N+1)])
        elif igrid == "r":
            sc = np.asarray([(k - N - .5)/N for k in range(1, N+1)])        
        # Get stretching function
        C = self.stretching(sc, Vstretching, theta_s, theta_b)        
        # Get depth of sigma layers
        z = self.get_zlev(H, C, hc, sc, zeta, Vtransform).transpose((2, 0, 1))        
        return z

    def stretching(self, sc, Vstretching, theta_s, theta_b):
        import numpy as np
        if Vstretching == 1:
            # Song and Haidvogel, 1994
            cff1 = 1.  / np.sinh(theta_s)
            cff2 = .5 / np.tanh(.5 * theta_s)
            C = (1.-theta_b) * cff1 * np.sinh(theta_s * sc) + \
                theta_b * (cff2 * np.tanh( theta_s * (sc + .5) ) - .5)
            return C
        
        if Vstretching == 2:
            # A. Shchepetkin (UCLA-ROMS, 2005) vertical stretching function        
            cff = 1 - sc**2
            Csur = ( 1. - np.cosh( theta_s * sc ) ) / ( np.cosh( theta_s ) - 1. )
            Cbot = -1. + ( np.sinh( theta_b * ( sc + 1. ) ) ) / ( np.sinh( theta_b ) )
            return cff * Csur + ( 1 - cff ) * Cbot        
        
        if Vstretching == 3:
            # R. Geyer BBL vertical stretching function.
            cff = .5 * ( 1. - np.tanh( 3. * ( sc + .5 ) ) )    
            Csur = -np.log( np.cosh( 3. * abs( sc )**theta_s)) / np.log( np.cosh( 3. ) )
            Cbot = np.log( np.cosh( 3. * ( sc + 1. )**theta_b)) / np.log( np.cosh( 3. ) ) - 1.
            return cff * Cbot + ( 1  - cff ) * Csur
    
        if Vstretching == 4:
            # A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
            if theta_s > 0:
                Csur = ( 1. - np.cosh(theta_s*sc) ) / ( np.cosh(theta_s) - 1. )
            else:
                Csur = -sc**2
    
            if theta_b > 0:
                Cbot = ( np.exp(theta_b*Csur)- 1. ) / ( 1. - np.exp(-theta_b) )
                return Cbot
            else:
                return Csur

    def get_zlev(self, h, sigma, hc, sc, ssh=0., Vtransform=2):
        import numpy as np
        # Get dimensions
        L, M = h.shape; N = len(sc)
        # Tiling bathymetry
        h = np.tile(h, (N, 1, 1)).transpose([1, 2, 0])
        # Tiling free-surface
        ssh = np.tile(ssh, (N, 1, 1)).transpose([1, 2, 0]) 
        # Tiling sigma-coordinate
        sc = np.tile(sc, (L, M, 1));   
        # Tiling stretching function
        sigma = np.tile(sigma, (L, M, 1));
        
        if Vtransform == 1: # ROMS 1999
            hinv = 1./h
            cff = hc * (sc - sigma)
            z0 = cff + sigma * h
            return z0 + ssh * (1. + z0*hinv)
        elif Vtransform == 2: # ROMS 2005
            z0 = ( hc * sc + sigma * h ) / ( h + hc )
            return ssh + (ssh + h) * z0
        
    def u2rho(self, u):
       import numpy as np
       """ interpolates from u-points to rho points """
       if u.ndim == 2:
          (ny,nx) = u.shape
          ur = np.ma.zeros((ny,nx+1),order = 'F')
          ur.mask = True
          ur[:,1:nx] = 0.5*(u[:,0:nx-1] + u[:,1:nx]);
          ur[:, 0] = ur[:, 1];
          ur[:,-1] = ur[:,-2];
       elif u.ndim == 3:
          (nz,ny,nx) = u.shape
          ur = np.ma.zeros((nz,ny,nx+1),order = 'F')
          ur.mask = True
          ur[:,:,1:nx] = 0.5*(u[:,:,0:nx-1] + u[:,:,1:nx]);
          ur[:,:, 0] = ur[:,:,1];
          ur[:,:,-1] = ur[:,:,-2];
       else:
          ur = 0
          print("error in u2rho")
    
       return ur     
    
    def v2rho(self, v):
       import numpy as np
       """ interpolates from v-points to rho points """    
       if v.ndim == 2:
          (ny,nx) = v.shape
          vr = np.ma.zeros((ny+1,nx),order = 'F')
          vr.mask = True
          vr[1:ny,:] = 0.5*(v[0:ny-1,:] + v[1:ny,:]);
          vr[ 0,:] = vr[1,:];
          vr[-1,:] = vr[-2,:];
       elif v.ndim == 3:
          (nz,ny,nx) = v.shape
          vr = np.ma.zeros((nz,ny+1,nx),order = 'F')
          vr.mask = True
          vr[:,1:ny,:] = 0.5*(v[:,0:ny-1,:] + v[:,1:ny,:]);
          vr[:, 0,:] = vr[:, 1,:];
          vr[:,-1,:] = vr[:,-2,:];
       else:
          vr = 0
          print("error in v2rho")
    
       return vr        
        
    def makecdf(self, cdf, bry, L, M, N, x, y, z, t0, t1, offset, var):  
        import numpy as np 
        """ Build output NetCDF """         
        f = netCDF4.Dataset(cdf, "w", format="NETCDF4")
        
        """ Create global attributes """
        f.creation_date = datetime.now().strftime("%d-%b-%Y %H:%M")
        f.northern_boundary = bry[0]
        f.southern_boundary = bry[1]
        f.eastern_boundary  = bry[2]
        f.western_boundary  = bry[3]
        
        """ Create dimensions """
        f.createDimension("y", L)
        f.createDimension("x", M)
        f.createDimension("z", N)
        f.createDimension("T", len(t0))
        f.createDimension("two", 2)
        f.createDimension("complex", 2)
        
        """ X """
        lon_rho = f.createVariable("lon_rho", "f4", \
                                   dimensions=("y", "x"), zlib=True)
        lon_rho.long_name = "longitude"
        lon_rho.units = "degree_east"
        lon_rho[:] = x
        
        """ Y """
        lat_rho = f.createVariable("lat_rho", "f4", \
                                   dimensions=("y", "x"), zlib=True)
        lat_rho.long_name = "latitude"
        lat_rho.units = "degree_north"
        lat_rho[:] = y
        
        """ Z """
        depth = f.createVariable("depth", "f4", \
                                 dimensions=("z"), zlib=True)
        depth.long_name = "depth of z-levels"
        depth.units = "meter"
        depth.positive = "up"
        if z:
            depth[:] = z
        
        """ T """
        time = f.createVariable("time", "f4", \
                                dimensions=("two", "T"), zlib=True)
        time.long_name = "user-selected periods"
        time.units = "seconds"
        time.offset = offset.strftime("%Y-%b-%d %H:%M")        
        time.first_column = "starting time"        
        time.second_column = "end time"
        time[:] = np.array([[(item - offset).total_seconds() for item in t0], \
                               [(item - offset).total_seconds() for item in t1]])
    
        """ Load metadata """
        # Keyword
        key = list(var)
        # Dimensions
        dim = [item[0] for item in var.values()]
        # Descriptions
        long = [item[2] for item in var.values()]
        # Units
        unit = [item[3] for item in var.values()]
        # Type
        tipo = [item[10] for item in var.values()]
        
        """ Loop along user-selected variables """
        for i1, i2, i3, i4, i5 in zip(key, dim, long, unit, tipo): 
            if i5 == "real":              
                A = f.createVariable(i1, "f4", dimensions=i2, \
                                     zlib=True)            
            else:
                A = f.createVariable(i1, "f4", dimensions=i2 + ("complex",), \
                                     zlib=True)
            A.long_name = i3
            if i4:
                A.units = i4
                
        f.close()
        
    def cdfseries(self, cdf, bry, N, T, z, TIME, offset, var):  
        import numpy as np
        """ Build output NetCDF """        
        f = netCDF4.Dataset(cdf, "w", format="NETCDF4")
        
        """ Create global attributes """
        f.creation_date = datetime.now().strftime("%d-%b-%Y %H:%M")
        f.northern_boundary = bry[0]
        f.southern_boundary = bry[1]
        f.eastern_boundary  = bry[2]
        f.western_boundary  = bry[3]
        
        """ Create dimensions """
        f.createDimension("z", N)
        f.createDimension("T", T) 
        f.createDimension("complex", 2)        
        
        """ Z """
        depth = f.createVariable("depth", "f4", dimensions=("z"), \
                                 zlib=True)
        depth.long_name = "depth of z-levels"
        depth.units = "meter"
        depth.positive = "up"
        if z:
            depth[:] = z
        
        """ T """
        time = f.createVariable("time", "f4", dimensions=("T"), \
                                zlib=True)
        time.long_name = "time"
        time.units = "seconds"
        time.offset = offset.strftime("%Y-%b-%d %H:%M")                     
        if TIME:
            time[:] = [(item - offset).total_seconds() for item in TIME]            
    
        """ Load metadata """
        # Keyword
        key = list(var)
        # Dimensions
        dim = [item[0] for item in var.values()]
        # Descriptions
        long = [item[2] for item in var.values()]
        # Units
        unit = [item[3] for item in var.values()]
        # Type
        tipo = [item[10] for item in var.values()]
        
        """ Loop along user-selected variables """
        for i1, i2, i3, i4, i5 in zip(key, dim, long, unit, tipo):
            if i5 == "real":              
                A = f.createVariable(i1, "f4", dimensions=i2[2:], \
                                     zlib=True)            
            else:
                A = f.createVariable(i1, "f4", dimensions=i2[2:] + ("complex",), \
                                     zlib=True)                     
            A.long_name = i3
            if i4:
                A.units = i4
            A[:] = np.nan
                
        f.close()
    
    def get_date(self, imn, emn):
        """ Get first and last date matching month serial numbers IMN and EMN """
        idate = self.dates[self.mon.index(imn)]
        edate = self.dates[len(self.mon) - 1 - self.mon[::-1].index(emn)]
        return idate, edate
    
    def util2(self, a, b, modo):
        import numpy as np
        import itertools
        
        if a.shape != b.shape:
            raise ValueError("Arrays must be the same size!")
        (L, M, T) = a.shape
        
        if modo == "max":
            ind = [item for sublist in np.argmax(a, axis=2).tolist() for item in sublist]
        elif modo == "min":
            ind = [item for sublist in np.argmin(a, axis=2).tolist() for item in sublist]
        else:
            raise ValueError("Mode must be either MAX or MIN")
        
        o = [x for x in itertools.product(list(range(L)), list(range(M)))]
        
        iidx = [x[0] for x in o]
        jidx = [x[1] for x in o]
    
        c = b[iidx, jidx, ind].reshape((L, M))
        return c
    
    def util3(self, a, b, modo):
        import numpy as np
        import itertools
        
        if a.shape != b.shape:
            raise ValueError("Arrays must be the same size!")
        (L, M, N, T) = a.shape
        
        if modo == "max":
            ind = [item for sublist in np.argmax(a, axis=3).tolist() for item in sublist]
            ind = [item for sublist in ind for item in sublist]
        elif modo == "min":
            ind = [item for sublist in np.argmin(a, axis=3).tolist() for item in sublist]
            ind = [item for sublist in ind for item in sublist]
        else:
            raise ValueError("Mode must be either MAX or MIN")
        
        o = [x for x in itertools.product(list(range(L)), list(range(M)), list(range(N)))]
        
        iidx = [x[0] for x in o]
        jidx = [x[1] for x in o]
        kidx = [x[2] for x in o]
    
        c = b[iidx, jidx, kidx, ind].reshape((L, M, N))
        return c
    
    def ProcessPeriod(self, v0, v1, n, mode):
        import numpy as np
        idx = np.ones(v1.shape, dtype=bool)
        if isinstance(v0, list):
            v1 = v1
            return v1, idx
        else:
            if mode == "avg":
                o = ( v0 * (n - 1) + v1 ) / n            
            elif mode == "min":
                o =  np.minimum(v0, v1)                
            elif mode == "max":
                o = np.maximum(v0, v1)                
        idx[o == v0] = False
        o = o
        return o, idx
        
if __name__ == "__main__" :  
    import sys
    orig = sys.stdout
    sys.stdout = open("log.txt", "w")
    root = Root()
    root.mainloop() 
    sys.stdout = orig