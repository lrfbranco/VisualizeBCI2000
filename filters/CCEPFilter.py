# importing various libraries
import numpy as np
from PyQt5 import QtWidgets
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg
from pyqtgraph.dockarea import *
import pyqtgraph.parametertree as ptree
from filters.filterBase.GridFilter import GridFilter
from base.SharedVisualization import saveFigure
from enum import Enum
from scipy.signal import find_peaks
from math import ceil
from nested_defaultdict_store import add_chunk, get_group, get_partial, to_py_type
from nested_defaultdict_store import root, metadata_keys
from pyqtgraph.parametertree import Parameter
from collections import deque # used for rolling buffer
import matplotlib.pyplot as plt # for plot testing
from scipy.signal import detrend

backgroundColor = (14, 14, 16)
highlightColor = (60, 60, 40)
highZValue = 1000

def parse_param(value):
      if value == "All":
          return None

      try:
          # if value is something like '130 (3)' or '1.5 (10)', split and parse first part
          numeric_part = value.split()[0]
          if '.' in numeric_part:
              return float(numeric_part)
          else:
              return int(numeric_part)
      except Exception as e:
          print(f"parse_param error: {e} | value: {value}")
          return value  # fallback, return original if parsing fails

class Column(Enum):
  Name = 0
  Electrode = 1
  Sig = 2
  AUC = 3
  RMS = 4
  P2P = 5

#taken from pyqtgraph example
## test add/remove
## this group includes a menu allowing the user to add new parameters into its child list
class ScalableGroup(ptree.parameterTypes.GroupParameter):
    def __init__(self, p, **opts):
        opts['type'] = 'group'
        opts['addText'] = "Add"
        opts['addList'] = ['str', 'float', 'int']
        ptree.parameterTypes.GroupParameter.__init__(self, **opts)
        self.p = p

        #enable/disable
        self.addChild({'name': 'Enable auto-detection', 'type': 'bool', 'value': 0})
        self.a = self.param('Enable auto-detection')
        self.a.sigValueChanged.connect(self.aChanged)
        #reference channel
        self.addChild({'name': 'Detection channel', 'type': 'str', 'value': "2", 'tip': 'Index or name'})
        self.b = self.param('Detection channel')
        self.b.sigValueChanged.connect(self.aChanged)
        self._dbsLayout = False # initialize DBS layout to unchecked
    def aChanged(self):
      self.p.setAutoDetect(self.a.value())
    def bChanged(self):
      self.p.setDetectChannels(self.b.value())
    
    def addNew(self, typ):
        val = {
            'str': '',
            'float': 0.0,
            'int': 0
        }[typ]
        self.addChild(dict(name="ScalableParam %d" % (len(self.childs)+1), type=typ, value=val, removable=True, renamable=True))

class TestBooleanParams(ptree.parameterTypes.GroupParameter):
  def __init__(self, p, **opts):
    self.p = p
    #opts['type'] = 'bool'
    #opts['value'] = True
    ptree.parameterTypes.GroupParameter.__init__(self, **opts)
    self.addChild({'name': 'Sort channels', 'type': 'bool', 'value': 0})
    self.a = self.param('Sort channels')
    self.a.sigValueChanged.connect(self.aChanged)

    self._dbsLeft = self.addChild({'name': 'DBS Layout Left', 'type': 'bool', 'value': 0})
    self._dbsLeft = self.param('DBS Layout Left')
    self._dbsLeft.sigValueChanged.connect(self.iChanged)

    self._dbsRight = self.addChild({'name': 'DBS Layout Right', 'type': 'bool', 'value': 0})
    self._dbsRight = self.param('DBS Layout Right')
    self._dbsRight.sigValueChanged.connect(self.jChanged)

    self.addChild({'name': 'Average CCEPS', 'type': 'bool', 'value': 0})
    self.b = self.param('Average CCEPS')
    self.b.sigValueChanged.connect(self.bChanged)

    self.addChild({'name': 'Threshold (STD)', 'type': 'slider', 'value': 2, 'span': np.linspace(0, 10, 100)})
    self.c = self.param('Threshold (STD)')
    #self.c.sigValueChanged.connect(self.cChanged)
    
    self.addChild({'name': 'Max Windows', 'type': 'int', 'value': 16, 'limits': [0, 100]})
    self.d = self.param('Max Windows')
    self.d.sigValueChanged.connect(self.dChanged)
            
    self.addChild({'name': 'Save Figures on Refresh', 'type': 'bool', 'value': 0})
    self.f = self.param('Save Figures on Refresh')
    self.f.sigValueChanged.connect(self.fChanged)
    
    self.addChild({'name': 'Save Figures', 'type': 'action'})
    self.g = self.param('Save Figures')
    self.g.sigActivated.connect(self.gChanged)

    self.addChild({'name': 'Clear Figures', 'type': 'action'})
    self.h = self.param('Clear Figures')
    self.h.sigActivated.connect(self.hChanged)

    self.addChild({'name': 'View ERNA Dictionary', 'type': 'action'})
    self.p.viewDictButton = self.param('View ERNA Dictionary')
    self.p.viewDictButton.sigActivated.connect(self.viewDictClicked)

    # Filtering buttons/dropdowns
    param = Parameter.create(name='Filters', type='group', children=[
        {'name': 'Frequency', 'type': 'list', 'limits': ['All', 100, 130, 160], 'value': 'All'},
        {'name': 'Amplitude', 'type': 'list', 'limits': ['All', 0.5, 1.0, 1.5], 'value': 'All'},
        {'name': 'Stimulation Channel', 'type': 'list', 'limits': ['All'] + [str(i) for i in range(1,17)], 'value': 'All'},
        {'name': 'Filter Data', 'type': 'action'},
        {'name': 'Plot Average ERNA', 'type': 'action'}
    ])
    self.addChild(param)

    self.filterFreq = self.param('Filters', 'Frequency')
    self.filterAmp = self.param('Filters', 'Amplitude')
    self.filterStim = self.param('Filters', 'Stimulation Channel')
    self.filterButton = self.param('Filters', 'Filter Data')
    self.plotAvg = self.param('Filters', 'Plot Average ERNA')

    self.filterButton.sigActivated.connect(self.local_filter_data)  # connect button to filter function
    self.plotAvg.sigActivated.connect(self.plot_avg_clicked)

  # maybe these can be removed for redundancy with local_filter_data ???
  def getFilterFreq(self):
    return self.filterFreq
  def getFilterAmp(self):
    return self.filterAmp
  def getFilterStim(self):
    return self.filterStim
  
  def viewDictClicked(self):
    freq = parse_param(self.filterFreq.value())
    amp = parse_param(self.filterAmp.value())
    stim = parse_param(self.filterStim.value())
    self.p.view_ERNA_dict(freq_filter=freq, amp_filter=amp, stim_filter=stim)

  def local_filter_data(self):
    freq = parse_param(self.filterFreq.value())
    amp = parse_param(self.filterAmp.value())
    stim = parse_param(self.filterStim.value())
    self.p.filter_data(freq, amp, stim)
  
  def plot_avg_clicked(self):
    freq = parse_param(self.filterFreq.value())
    amp = parse_param(self.filterAmp.value())
    stim = parse_param(self.filterStim.value())
    self.p.plot_average_erna(freq_filter=freq, amp_filter=amp, stim_filter=stim)

  def aChanged(self):
    self.p.setSortChs(self.a.value())
  def bChanged(self):
    self.p.setAvgPlots(self.b.value())
  def cChanged(self):
    self.p.setStdDevState(self.c.value())
  def dChanged(self):
    self.p.setMaxWindows(self.d)
  def eChanged(self):
    self.p.setMaxPlots(self.e)
  def fChanged(self):
    self.p.setSaveFigs(self.f.value())
  def gChanged(self):
    self.p.saveFigures()
  def hChanged(self):
    self.p.clearFigures()
  def iChanged(self):
    # avoid sorting channels & dbs simultaneously
    if self.p.setSortChs:
        self.p.setSortChs = False
        #QtWidgets.QMessageBox.warning(None, "Warning", "Sort Channels has been disabled because DBS Layout was activated.")
    self.p.applyDBSLayout(self._dbsLeft.value() or self._dbsRight.value())
  def jChanged(self):
    if self.p.setSortChs:
        self.p.setSortChs = False
        #QtWidgets.QMessageBox.warning(None, "Warning", "Sort Channels has been disabled because DBS Layout was activated.")
    self.p.applyDBSLayout(self._dbsLeft.value() or self._dbsRight.value())


class CCEPFilter(GridFilter):
  def __init__(self, area, bciPath, stream):
    super().__init__(area, bciPath, stream)
    self.aucThresh = 0
    self.numTrigs = 0

    self.erna_buffer = deque(maxlen=3) # create a rolling buffer of the last 10 erna detections

  def publish(self):
    super().publish()

    self.pens = [pg.mkPen(x) for x in np.linspace(0, 1, 256)] #create all the pens we could ever need
    self.gridNums = pg.GraphicsLayoutWidget(title="CCEP Aggregate")
    self.table = QtWidgets.QTableWidget()

    settingsD = Dock("Settings")
    settingsLab = QtWidgets.QLabel("Settings", objectName="h1")

    #create parameter tree
    params = [
        TestBooleanParams(name= 'General Options', p=self, showTop=False),
        ScalableGroup(name="Auto Detect Options", p=self, tip='Click to add channels'),
    ]

    self.p = ptree.Parameter.create(name="Settings", type='group', children=params, title=None)
    self.t = ptree.ParameterTree()
    self.t.setParameters(self.p)

    self.epoch_count = 0
    self.p.addChild({'name': 'Epoch Count', 'type': 'int', 'value': self.epoch_count, 'readonly': True})
    self.epochDisplay = self.p.param('Epoch Count')
    self.epochDisplay.setValue(self.epoch_count)

    #self.t.setParameters(params)
    #self.t.resizeColumnToContents(0)
    self.t.header().setSectionResizeMode(pg.QtWidgets.QHeaderView.Stretch)
    #self.t.show()
    settingsD.addWidget(settingsLab)
    settingsD.addWidget(self.t)

    d2 = Dock("Table", widget=self.table)
    self.area.addDock(settingsD)
    self.area.addDock(d2, position='above', relativeTo=settingsD)
    self.area.addDock(Dock("Plots", widget=self.gridPlots), position='above', relativeTo=d2)

  def loadSettings(self):
    super().loadSettings()
    self._stds = self.settings.value("stds", 5)

    self._maskStart = self.settings.value("maskStart", -5)
    self._maskEnd = self.settings.value("maskEnd", 15)

    pState = self.settings.value("pState", {})
    # if hasattr(pState, 'children'):
    #   children = pState['children']

    if pState != {}:
      self.p.restoreState(pState, addChildren=False, removeChildren=False)
    self._maxWindows = self.p.child('General Options')['Max Windows']
    self._sortChs = self.p.child('General Options')['Sort channels']
    self._trigCh = self.p.child('Auto Detect Options')['Detection channel']
    self.epoch_count = 0
    self.epochDisplay.setValue(self.epoch_count)    # set epoch count to zero upon restarting GUI

    # reset filter dropdowns to be empty upon restarting GUI
    # test_params = self.p.child('General Options')
    # freq_param = test_params.getFilterFreq()
    # amp_param = test_params.getFilterAmp()
    # stim_param = test_params.getFilterStim()
    # freq_defaults = ['All']
    # amp_defaults = ['All']
    # stim_defaults = ['All']
    # freq_param.setLimits(freq_defaults)
    # amp_param.setLimits(amp_defaults)
    # stim_param.setLimits(stim_defaults)
    # freq_param.setValue('All')
    # amp_param.setValue('All')
    # stim_param.setValue('All')

  def saveSettings(self):
    super().saveSettings()

    self.settings.setValue("maskStart", self._maskStart)
    self.settings.setValue("maskEnd", self._maskEnd) 
    self.settings.setValue("pState", self.p.saveState())

  #an attempt to abstract plotting from BCI2000
  def checkPlot(self):
    return self.comm.evaluate("CCEPTriggered")
  
  #define shared states we need for filter
  @property
  def sharedStates(self):
    return ["CCEPTriggered", "StimulatingChannel"]
  
  #define abstract methods
  def receiveStates(self, state):
    #get CCEPTriggered state to detect CCEPs
    self.numTrigs += np.count_nonzero(state[0]) # keep running count of when CCEPTrigger occurs
    # print(f"[DEBUG] numTrigs: {self.numTrigs}")   # confirm numTrigs is updating

    #find stim ch if possible
    if np.shape(state)[0] > 1:
      stimCh = state[1].nonzero()[0]
      if stimCh.size > 0:
        #just get first non-zero value
        chBits = state[1][stimCh[0]]
        self.stimChs.clear()
        chBinary = str("{0:b}".format(chBits))
        for b in range(len(chBinary)): #32 bit state
          if chBinary[len(chBinary) - b - 1] == '1':
            #print(self.chNames[b] + " at " + str(b))
            if b < len(self.chNames):
              self.stimChs.append(self.chNames[b]) #append ch name
            else:
              print(f"[WARN] Bit position {b} exceeds available channels ({len(self.chNames)})") #append ch name

  def plot_average_erna(self, freq_filter, amp_filter, stim_filter):
    """
    For each channel:
    1) pull out all stored epochs macthing parameters & channel name
    2) compute mean and SEM across those epochs
    3) plot in 4x4 layout (same as raw data layout)
    """
    results = self.filter_data(freq_filter, amp_filter, stim_filter)
    if not results:
      print("No matching epochs found.")
      return
    
    sample_count = results[0]['data'].shape[0]
    time_axis = np.linspace(-self.baselineLength, self.ccepLength, sample_count)  # time for each captured epoch

    avgPlotMap = {}

    avg_dock = Dock("Average ERNA", size=(1, 1), closable=True)
    avg_widget = pg.GraphicsLayoutWidget()
    avg_dock.addWidget(avg_widget)
    self.area.addDock(avg_dock, position='right', relativeTo=None)

    for chName, (row, col) in self._getChannelPositions().items():
      ch_epochs = [chunk['data'] 
                    for chunk in results 
                    if chunk.get('channel') == chName]   # use stored metadata channel name
      if not ch_epochs:
          continue

      data_array = np.stack(ch_epochs, axis=0)
      avg = data_array.mean(axis=0)
      sem = data_array.std(axis=0) / np.sqrt(data_array.shape[0])

      p = pg.PlotItem(title=f"{chName}\n(n={data_array.shape[0]})")
      p.plot(time_axis, avg, pen=pg.mkPen('b', width=1.5))
      upper = p.plot(time_axis, avg + sem, pen=None)
      lower = p.plot(time_axis, avg - sem, pen=None)
      fill = pg.FillBetweenItem(upper, lower, brush=(0, 0, 255, 50))
      p.addItem(fill)

      p.setTitle(f"{chName}\n(n={data_array.shape[0]})")
      p.setLabel('bottom', 'Time (ms)')
      p.setLabel('left', 'Amplitude (µV)')
      p.setYRange(-1000, 1000)

      avgPlotMap[chName] = p

    self._applyLayout(avg_widget, avgPlotMap)
    avg_widget.show()

  def plot(self, data):
    #if self.checkPlot():
    if self.numTrigs > 0:
      # print("plotting")
      self.numTrigs -= 1

      #process data
      aocs = []
      chunk = False
      peaks = None
      if self.p.child('Auto Detect Options')['Detection channel']:
        #trigCh = self.p.child('General Options')['Sort channels']
        #get channel to use as trigger
        try:
          chIndex = int(self._trigCh)
        except ValueError:
          try:
            #not index, gotta be ch name
            chIndex = self.chNames.index(self._trigCh)
          except:
            self.logPrint(self._trigCh + " is not a valid channel name")
            chIndex = 1
        self.trigData = data[chIndex - 1] #convert to 0-based, retreive last channel (trigger channel)

        #get chunks by peaks
        #hard code parameters just to test
        # print(f"Trigger data max: {np.max(self.trigData)}, min: {np.min(self.trigData)}")

        peaks, properties = find_peaks(self.trigData, 1, distance=1)    # changed values for testing
        # print(self.trigData)
        # print(f"Found {len(peaks)} peaks") # TEMPORARILY DISABLE FOR TESTING
        if len(peaks) >= 1:   # our data has only one peak, so consider it
          chunk = True
        else:
          chunk = False
        # print(chunk)

        # #chunk based off artifact
        # for i, ch in enumerate(self.chTable.values()):
        #   ch.chunkData(data[i], peaks) #chunks and computes
        #   aocs.append(ch.auc)
      
      #compute and chunk data
      avgPlots = self.p.child('General Options')['Average CCEPS']
      if chunk:
        # Generate metadata ONCE per stimulation event
        fake_meta = {
            'amplitude': np.random.choice([0.5, 1.0, 1.5]),
            'frequency': np.random.choice([100, 130, 160]),
            'stim_channel': np.random.choice(list(range(1,17))),
            'trial_id': self.epoch_count
        }

        # process all channels
        for i, ch in enumerate(self.chTable.values()):
            ch.chunkData(data[i], peaks, avgPlots)
            last_peak_idx = peaks[-1]
            last_peak_time = last_peak_idx / self.sr * 1000 # in ms
            ch.last_peak_time = last_peak_time
            ch.computeFeatures()
            meta = fake_meta.copy()
            meta['channel'] = ch.title if hasattr(ch, 'title') else self.chNames[i]
            add_chunk(ch.data.copy(), meta)

            is_detected = ch.auc > 1 # arbitrary threshold used for testing
            self.erna_buffer.append({
              'detected': is_detected,
              'auc': ch.auc
            })

        self.epoch_count += 1
        self.epochDisplay.setValue(self.epoch_count)
        self.update_filter_dropdowns()
        # stability = self.compute_stability()
        # if stability:
        #     print(f"ERNA Stability ➔ Detection rate: {stability['detection_rate']*100:.1f}% | Mean AUC: {stability['mean_auc']:.2f} | Variance: {stability['var_auc']:.2f}")


      else:
          for i, ch in enumerate(self.chTable.values()):
              ch.computeData(data[i], avgPlots)  # compute data
              last_peak_idx = peaks[-1]
              last_peak_time = last_peak_idx / self.sr * 1000 # in ms
              ch.last_peak_time = last_peak_time
              ch.computeFeatures()
              # print("we are in the else")
      for i, ch in enumerate(self.chTable.values()):
          aocs.append(ch.auc)      
      #send processed data
      self.dataProcessedSignal.emit(aocs)

      #we scale by 10 cause slider can only do ints
      stds = self.p.child('General Options')['Threshold (STD)']
      self.aucThresh = np.std(aocs) * stds

      #plot!
      self._renderPlots()

  def setConfig(self):
    super().setConfig()

    self.gridNums.clear()

    self.chPlot = list(range(self.channels))
    self.tableRows = list(range(self.channels))
    #self.chPlot = {}
    self.chTable = {}
    self.regs = list(range(self.channels))
    #init variables
    self.baselineLength = self.getParameterValue("BaselineEpochLength")
    self.latStart = 0
    self.latStartSamples = self._maskStart
    self.ccepLength = self.getParameterValue("CCEPEpochLength")
    self.sr = self.getParameterValue("SamplingRate")
    self.baseSamples = self.msToSamples(self.baselineLength)
    self.ccepSamples = self.msToSamples(self.ccepLength)
    self.trigSamples = self._maskEnd 
    self.trigLatLength = self.trigSamples * 1000.0 / self.sr
    
    #redefine element size
    self.elements = self.baseSamples + self.ccepSamples
    self.x = np.linspace(-self.baselineLength, self.ccepLength, self.elements)

    #to visualize stimulating channels if we can
    self.stimChs = []

    #go thru all channels for table
    count = 0
    for chName in self.chNames:
      sub1 = self.gridNums.addLayout()
      sub1.addLabel("<b>%s"%(chName), size='20pt', bold=True)

      sub1.nextRow()
      self.chTable[chName] = CCEPCalc(self, ch=count, title=chName)
      count = count + 1

    #only initialize plots up to max number 
    for r in range(self.numRows):
      for c in range(self.numColumns):
        ch = r*self.numColumns+c
        #print(self.chPlot)
        if ch < self.windows:
          chName = self.chNames[ch]
          self.chPlot[ch] = CCEPPlot(self, title=chName, row=self.chTable[chName])
          self.gridPlots.addItem(self.chPlot[ch])
          if ch != 0:
            self.chPlot[ch].setLink(self.chPlot[ch-1])
          else:
            self.chPlot[ch].showAxis('left')
            self.chPlot[ch].showAxis('bottom')
      self.gridPlots.nextRow()
    
    #give first plot a friend
    if self.windows > 1:
      self.chPlot[0].friend = self.chPlot[self.windows-1] #give first plot a friend
    
    #table
    self.table.setRowCount(self.channels)
    self.table.setColumnCount(len(Column))
    self.table.setHorizontalHeaderLabels([c.name for c in Column])
    for ch, chName in enumerate(self.chNames):
      c = ch % self.numColumns
      r = int(np.floor(ch/self.numColumns))
      #add table
      if self.elecDict:
        eName = self.elecDict[chName]
      else:
        eName = ""
      self.tableRows[ch] = TableRow(ch, r, c, self.chPlot[ch], chName, eName, self.channels)
      self.tableRows[ch].addRow(self.table)
      
      if ch < self.windows:
        self.chPlot[ch].setRow(self.tableRows[ch]) #link plot to row
    
    #finish table
    if not self.elecDict:
      self.table.setColumnHidden(Column.Electrode.value, True) #hide electrode name col
    else:
      self._hideNonElectrodes()
    self.table.resizeColumnsToContents()
    #self.table.setSortMode(0,'index')
    self.table.sortItems(0, QtCore.Qt.DescendingOrder)
    for i, chName in enumerate(self.chNames):
      self.chTable[chName].setTableItem(i)
    #self.setSortChs(self._sortChs)
    #make sure user can't change sorting
    for i in range(self.table.columnCount()):
      h = self.table.horizontalHeaderItem(i)
      #h.setFlags(h.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
      self.table.setHorizontalHeaderItem(i, h)
    # print(self.table.horizontalHeaderItem(1).text())
    self.selectedRows = []
    #self.table.itemClicked.connect(self.tableItemClickedCallback)
    self.table.itemSelectionChanged.connect(self._selectionChanged)
    #self.table.itemChanged.connect(self._itemChanged)

  def setStdDevState(self, value):
    self._stds = value
    #self.stdSpin.setToolTip(str(value/10))
  def setMaxWindows(self, spin):
    self._maxWindows = spin.value()
  def setAvgPlots(self, state):
    self._avgPlots = state

  
  def setSortChs(self, state):
    if hasattr(self, 'chTable') and self._sortChs and not state:
      #re-initialize order
      for ch in self.chTable.values():
        ch.totalChanged(False)
      self.table.sortItems(Column.Name.value, QtCore.Qt.DescendingOrder)
    self._sortChs = state
    if hasattr(self, 'chTable'):
        self._renderPlots(newData=False)

  def _getChannelPositions(self):
      """
      Returns a dict {channelName: (row, col)} based on current DBS or default layout.
      """
      leftOn  = self.p.param('General Options', 'DBS Layout Left').value()
      rightOn = self.p.param('General Options', 'DBS Layout Right').value()

      if leftOn or rightOn:
          # Define left and right DBS mappings
          left_positions = {
              "Lv4_L":  (0, 1), "Lv3C_L": (1, 0), "Lv3A_L": (1, 1), "Lv3B_L": (1, 2),
              "Lv2C_L": (2, 0), "Lv2A_L": (2, 1), "Lv2B_L": (2, 2), "Lv1_L":  (3, 1),
          }
          right_positions = {
              "Lv4_R":  (0, 5), "Lv3C_R": (1, 4), "Lv3A_R": (1, 5), "Lv3B_R": (1, 6),
              "Lv2C_R": (2, 4), "Lv2A_R": (2, 5), "Lv2B_R": (2, 6), "Lv1_R":  (3, 5),
          }
          positions = {}
          if leftOn:
              positions.update(left_positions)
          if rightOn:
              positions.update(right_positions)
          return positions
      else:
          # Default row-major layout
          positions = {}
          for i, ch in enumerate(self.chNames):
              row = i // self.numColumns
              col = i % self.numColumns
              positions[ch] = (row, col)
          return positions

  def _applyLayout(self, layoutWidget, plotMap):
      """
      Clear existing items and place each PlotItem from plotMap into layoutWidget
      using positions from _getChannelPositions().
      """
      # Remove all existing items
      layoutWidget.clear()
      # Fetch where each channel should go
      positions = self._getChannelPositions()
      # Add only valid PlotItems
      for chName, (r, c) in positions.items():
          item = plotMap.get(chName)
          if isinstance(item, pg.PlotItem):
              layoutWidget.addItem(item, row=r, col=c)

  def applyDBSLayout(self, state):
      """
      Apply DBS or default layout to the raw-data gridPlots.
      """
      if not hasattr(self, 'chTable'):
          return

      # 1) Update table to show only channels in active layout
      self.updateTableForDBSLayout()

      # 2) Build mapping: channelName -> PlotItem (filter out any non-PlotItems)
      plotMap = {}
      for idx, chName in enumerate(self.chNames):
          item = self.chPlot[idx]
          if isinstance(item, pg.PlotItem):
              plotMap[chName] = item

      # 3) Re-layout the gridPlots area
      self._applyLayout(self.gridPlots, plotMap)

      # 4) Update flags and re-render plots
      self._dbsLayout = state
      self._renderPlots(newData=False)

  def updateTableForDBSLayout(self):
    leftOn = self.p.param('General Options', 'DBS Layout Left').value()
    rightOn = self.p.param('General Options', 'DBS Layout Right').value()

    if not leftOn and not rightOn:
        for row in range(self.table.rowCount()):
            self.table.setRowHidden(row, False)
        return

    for row in range(self.table.rowCount()):
      ch_item = self.table.item(row, Column.Name.value)
      if not ch_item:
        continue
      ch_name = ch_item.text()
      is_left = ch_name.endswith('_L')
      is_right = ch_name.endswith('_R')

      visible = (
        (leftOn and is_left) or (rightOn and is_right)
      )
      self.table.setRowHidden(row, not visible)

  def compute_stability(self):
    if not self.erna_buffer:
      return None
    
    detections = [entry['detected'] for entry in self.erna_buffer]
    aucs = [entry['auc'] for entry in self.erna_buffer]

    detection_rate = sum(detections) / len(detections) # number of times erna was detected over the last 10 stimulations
    mean_auc = np.mean(aucs)
    var_auc = np.var(aucs)

    return {
        'detection_rate': detection_rate,
        'mean_auc': mean_auc,
        'var_auc': var_auc
    }

  def filter_data(self, freq, amp, stim):
      # create query dictionary
      query = {}
      if freq != 'All':
          query['frequency'] = freq
      if amp != 'All':
          query['amplitude'] = amp
      if stim != 'All':
          query['stim_channel'] = stim
      
      print(f"\nFilter query → {query}")

      # retrieve matching data
      results = get_partial(query)

      if results:
          print(f"Number of Matching Results = {len(results)}")
      else:
          print("No matching results found.")
      return results
      #print(f"Filter results: {len(results)} entries found for query {query}")
  
  def update_filter_dropdowns(self):
    test_params = self.p.child('General Options')

    freq_param = test_params.getFilterFreq()
    amp_param = test_params.getFilterAmp()
    stim_param = test_params.getFilterStim()

    current_freq = parse_param(freq_param.value())
    current_amp = parse_param(amp_param.value())
    current_stim = parse_param(stim_param.value())

    freq_counts = {}
    amp_counts = {}
    stim_counts = {}

    for freq in [100, 130, 160]:
        results = get_partial({'frequency': freq})
        if results:
            freq_counts[freq] = len(set(entry['trial_id'] for entry in results))

    for amp in [0.5, 1.0, 1.5]:
        results = get_partial({'amplitude': amp})
        if results:
            amp_counts[amp] = len(set(entry['trial_id'] for entry in results))

    for stim in range(1, 17):
        results = get_partial({'stim_channel': stim})
        if results:
            stim_counts[stim] = len(set(entry['trial_id'] for entry in results))

    freq_options = ['All'] + [f"{freq} ({count})" for freq, count in freq_counts.items()]
    amp_options = ['All'] + [f"{amp} ({count})" for amp, count in amp_counts.items()]
    stim_options = ['All'] + [f"{stim} ({count})" for stim, count in stim_counts.items()]

    freq_param.setLimits(freq_options)
    amp_param.setLimits(amp_options)
    stim_param.setLimits(stim_options)

    def restore_selection(param, current_val, options):
        for opt in options:
            if parse_param(opt) == current_val:
                param.setValue(opt)
                return
        param.setValue('All')

    restore_selection(freq_param, current_freq, freq_options)
    restore_selection(amp_param, current_amp, amp_options)
    restore_selection(stim_param, current_stim, stim_options)

  def setSaveFigs(self, state):
    self._saveFigs = state
  def setAutoDetect(self, state):
    self._autoDetect = state
    #self.trigChForm.setReadOnly(not state)
  def setDetectChannels(self, text):
    print(text)
    self._trigCh = text
      
  def msToSamples(self, lengthMs):
    return int(ceil(lengthMs * self.sr/1000.0))

  def updateParameter(self, latStart, newLat):
    if newLat != self.trigLatLength:
      self.trigLatLength = newLat
      self.trigSamples = self.msToSamples(newLat)
      self._maskEnd = self.latStartSamples + self.trigSamples
    if latStart != self.latStart:
      self.latStart = latStart
      self.latStartSamples = self.msToSamples(latStart)
      self._maskStart = self.latStartSamples
  
  def clearFigures(self):
    if self.p.child('General Options')['Save Figures on Refresh']:
      self.saveFigures()
    for i in range(0, self.windows):
      children = self.chPlot[i].listDataItems()
      for child in children[1:]: #save first plot
        self.chPlot[i].removeItem(child)
    children[0].setPen(pg.mkPen('b')) #blend in

    for t in self.chTable.values():
      t.totalChanged(0)
      t.database = []

  def saveFigures(self):
    #beautiful hack
    #mimic minimal mouse click
    class Dummy():
      def __init__(self):
        acceptedItem = None
      def screenPos(self):
        return pg.QtCore.QPointF(0,0)

    if self.gridPlots.getItem(0,0) != None:
      scene = self.gridPlots.scene()
      if not hasattr(scene, "contextMenuItem"):
        vb = self.gridPlots.getItem(0,0).getViewBox()
        
        event = Dummy()
        event.acceptedItem = vb
        vb.raiseContextMenu(event)
      scene.showExportDialog()

    return
    saveFigure(self.gridPlots)
    
  def _renderPlots(self, newData=True):
    #update table with new data
    for ch in self.chTable.values():
      ch.totalChanged(ch.auc > self.aucThresh)
    
    #sort table with updated numbers, if toggled
    if self._sortChs:
      self.table.sortItems(Column.Name.value, QtCore.Qt.DescendingOrder)
    
    #plot
    i = 0
    for r in range(self.table.rowCount()):
      if i == self.windows:
        break
      if not self.table.isRowHidden(r):
        chName = self.table.item(r, Column.Name.value).text()
        self.chPlot[i].changePlot(chName, self.chTable[chName])
        self.chPlot[i].plotData(newData)
        i+=1

  def _changeBackgroundColor(self, row, emph):
    if row >= self.windows:
      return
    c = backgroundColor
    if emph:
      c = highlightColor
    self.chPlot[row].vb.setBackgroundColor(c)
    chName = self.table.item(row,Column.Name.value).text()
    self.chPlot[row].selected = emph
    self.chTable[chName].selected = emph
  
  def _selectionChanged(self):
    items = self.table.selectedItems()
    newRows = []
    for p in items:
      if p.row() not in self.selectedRows:
        self._changeBackgroundColor(p.row(),True)
      newRows.append(p.row())
    for oldS in self.selectedRows:
      if oldS not in newRows:
        self._changeBackgroundColor(oldS, False)
    self.selectedRows = newRows

  ####---inherited slots---####
  def acceptElecNames(self, elecDict):
    super().acceptElecNames(elecDict)
    if hasattr(self, 'chTable'):
      for name in self.chNames:
        r = self.table.row(self.chTable[name].tableItem)
        self.table.item(r,Column.Electrode.value).setData(QtCore.Qt.DisplayRole, self.elecDict[name])
      self.table.setColumnHidden(Column.Electrode.value, False)
      self._hideNonElectrodes()
    else:
      #we haven't initialized table yet
      pass
  
  #if we have electrode names, we can hide channels that aren't electrodes
  def _hideNonElectrodes(self):
    if self.elecDict and hasattr(self, 'chTable'):
      for r in range(self.table.rowCount()):
        if self.table.item(r,Column.Electrode.value).text() == "":
          self.table.setRowHidden(r, True)
        else:
          self.table.setRowHidden(r, False)

  def view_ERNA_dict(self, freq_filter, amp_filter, stim_filter):
    print("\n==== Stored ERNA Dictionary ====")

    def traverse(node, freq=None, amp=None, stim=None):
        if '_chunks' in node:
            trial_ids = set()
            for entry in node['_chunks']:
                trial_id = entry.get('trial_id')
                trial_ids.add(trial_id)
            num_epochs = len(trial_ids)

            # apply filters: if any filter is set and doesn't match, skip
            if ((amp_filter is not None and amp != amp_filter)):
                return
            if (freq_filter is not None and freq != freq_filter):
                return
            if (stim_filter is not None and stim != stim_filter):
                return
            print(f"Amplitude: {amp} mA | Frequency: {freq} Hz | Stim Channel: {stim} | Number of epochs: {num_epochs}")
        else:
            for key, subnode in node.items():
                if amp is None:
                    traverse(subnode, freq=freq, amp=key, stim=stim)
                elif freq is None:
                    traverse(subnode, freq=key, amp=amp, stim=stim)
                elif stim is None:
                    traverse(subnode, freq=freq, amp=amp, stim=key)
                else:
                    traverse(subnode, freq=freq, amp=amp, stim=stim)

    traverse(root)

    print("================================\n")

class TableRow():
  class TableItem(QtWidgets.QTableWidgetItem):
    def __init__(self, p, parent=0):
      super().__init__(parent)
      self.p = p
      #make un-editable
      self.setFlags(self.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
      pass
    #define less than (<) operator for table sorting
    def __lt__(self, b):
      return bool((self.p.rank + self.p.sig*self.p.max) < (b.p.rank + b.p.sig*b.p.max))   # force comparison result to be native python bool to avoid type error
  def __init__(self, i, row, col, fig, chName, elName, maxVal):
    self.plotR = row
    self.plotC = col
    self.rank = maxVal - i
    self.sig = 0
    self.max = maxVal
    self.oldRow = i
    self.chName = chName
    self.elName = elName
    self.figure = fig
  def getRowNumber(self):
    return self.n.row()
  def addRow(self, table):
    self.n = self.TableItem(self, self.chName) #save order as rank
    s = self.TableItem(self)
    a = self.TableItem(self)
    e = self.TableItem(self, self.elName)
    r = self.TableItem(self)
    p = self.TableItem(self)

    table.setItem(self.oldRow, Column.Name.value, self.n)
    table.setItem(self.oldRow, Column.Electrode.value, e)
    table.setItem(self.oldRow, Column.Sig.value, s)
    table.setItem(self.oldRow, Column.AUC.value, a)
    table.setItem(self.oldRow, Column.RMS.value, r)
    table.setItem(self.oldRow, Column.P2P.value, p)

class CCEPPlot(pg.PlotItem):
  def __init__(self, parent, title, row):
    super().__init__(title=title)
    self.p = parent
    self.name = title
    self.link = row
    self.selected = False
    self.tableRow = None
    #prepare view
    axView = self.getViewBox()
    axView.disableAutoRange()
    axView.setMouseEnabled(x=False, y=True)
    axView.setDefaultPadding(0)
    xLim = -self.p.baselineLength
    yLim = self.p.ccepLength
    axView.setXRange(xLim, yLim, padding=0)
    axView.setYRange(-1000, 1000)

    self.setMinimumSize(100,100)

    #stim artifact filter
    self.latLow = self.p.latStartSamples*1000.0/self.p.sr
    self.latHigh = self.p.trigLatLength
    self.latReg = pg.LinearRegionItem(values=(self.latLow, self.latHigh), movable=True, brush=(9, 24, 80, 100), 
                                      pen=pg.mkPen(color=(9, 24, 80), width=1, style=QtCore.Qt.DotLine), bounds=[xLim, yLim])
    self.latReg.setZValue(highZValue) #make sure its in front of all plots
    #callbacks
    self.latReg.sigRegionChanged.connect(self.regionChanged)
    self.addItem(self.latReg)

    #initialize average plot
    self.avg = self.plot(x=self.p.x, y=np.zeros(self.p.elements), pen=pg.mkPen(backgroundColor)) #filler data
    self.avg.setZValue(highZValue-1) #behind filter region, in front of every other plot

    #self.backgroundC = (14, 14, 14)
    self.vb.setBackgroundColor(backgroundColor)
  def setRow(self, row):
    self.tableRow = row
  
  def setLink(self, plt):
    self.friend = plt #each plot gets one friend they affect
    self.getViewBox().setYLink(self.friend.getViewBox())

  def regionChanged(self, reg):
    newReg = reg.getRegion()
    self.latHigh = newReg[1]
    self.latLow = newReg[0]
    if self.latHigh != self.friend.latHigh or self.latLow != self.friend.latLow:
      self.p.updateParameter(self.latLow, self.latHigh) 
      self.friend.latReg.setRegion(reg.getRegion())
  
  def changePlot(self, name, link):
    #have we changed channels
    if self.name != name:
      self.setTitle(name)
      self.name = name
      self.link = link #update link
      if len(self.link.database) > 0:
        #change data of all plots but average
        for f, d in zip(self.listDataItems()[1:], self.link.database):
          f.setData(x=self.p.x, y=d, useCache=True)

      #change background based on selected
      if self.selected and not self.link.selected:
        self.vb.setBackgroundColor(backgroundColor)
        self.selected = False
      elif not self.selected and self.link.selected:
        self.vb.setBackgroundColor(highlightColor)
        self.selected = True

  #plot: new name and link is only considered if we are dynamically sorting
  def plotData(self, newData, maxPlots=0):
    if newData:
      children = self.listDataItems() #all plots
      
      if len(self.link.database) > len(children):
        #update plots with newly computed data (will be true if chunking)
        newPlots = len(self.link.database) - len(children)
        #print("new plots: " + str(newPlots))
        for i in range(newPlots):
          self.plot(x = self.p.x, y = self.link.database[-i - 1])

      expColors = [(255) * (1 - 2**(-x)) for x in np.linspace(0+1/(len(children)+1), 1-1/(len(children)+1), len(children))]
      for child, c in zip(children[1:], expColors):
        child.setPen(self.p.pens[int(c)])

      #plot new data
      p2 = 255*(1-2.5**(-1*(1-1/(len(self.link.database)+1))))
      self.plot(x=self.p.x, y=self.link.database[-1], useCache=True, pen=self.p.pens[int(p2)], _callSync='off')
    
    #update average plot
    if self.link.significant:
      p = pg.mkPen('y', width=1.5) #ccep!
    else:
      p = pg.mkPen('w', width=1.5)
    if self.name in self.p.stimChs:
      p = pg.mkPen('c', width=1.5)
    self.avg.setData(x=self.p.x, y=self.link.data, useCache=True, pen=p)

class CCEPCalc():
  def __init__(self, parent, ch, title):
    self.p = parent
    self.ch = ch

    self.significant = False
    self.selected = False
    self.database = []
    self.auc = 0
    self.rms = 0
    self.p2p = 0
    self.data = np.zeros(self.p.elements)

  def getActiveData(self, data):
    p1 = data[:self.p.baseSamples+self.p.latStartSamples]
    p2 = data[self.p.baseSamples+self.p.trigSamples:]
    d = np.concatenate((p1, p2))
    #print(np.shape(d))
    return p2

  def setTableItem(self, ch):
    self.tableItem = self.p.table.item(ch,Column.Name.value)

  #t = boolean, if significant or not
  def totalChanged(self, t):
    empColor = pg.QtGui.QColor(56,50,0)
    self.tableItem.p.sig = t
    r = self.p.table.row(self.tableItem) #find new row we are at
    self.p.table.item(r,Column.Sig.value).setData(QtCore.Qt.DisplayRole, int(t))
    self.p.table.item(r,Column.AUC.value).setData(QtCore.Qt.DisplayRole, int(self.auc))
    self.p.table.item(r,Column.RMS.value).setData(QtCore.Qt.DisplayRole, round(self.rms))
    self.p.table.item(r,Column.P2P.value).setData(QtCore.Qt.DisplayRole, int(self.p2p))
    self.significant = t

    if self.significant and self.tableItem.background() != empColor:
      self.tableItem.setBackground(empColor)
    elif not self.significant and self.tableItem.background() == empColor:
      self.tableItem.setBackground(pg.QtGui.QColor(0,0,0,0)) #transparent

  def computeData(self, newData, avgPlots=True):        
    #new data, normalize amplitude with baseline data
    if self.p.baseSamples == 0:
      self.data = newData
      #stdBase = 0
    else:
      #avBase = np.mean(newData[:self.p.baseSamples])
      avBase = np.median(newData[:self.p.baseSamples])
      #stdBase = np.std(self.rawData[:self.p.baseSamples], dtype=np.float64)
      self.data = np.subtract(newData, avBase)

    #store data
    self.database.append(self.data.copy())

    #possibly change to average, before we detect ccep
    if avgPlots:
      #calculate average of plots
      self.data = np.mean(self.database, axis=0)
    
    #get area under the curve
    ccepData = self.getActiveData(self.data)
    if ccepData.size > 0: # ensure array is not empty
      normData = ccepData - np.mean(ccepData)
      self.auc = np.trapz(abs(normData))/1e3
    else:
      self.auc = 0
    

  def chunkData(self, newData, peaks, avgPlots=True):
    for peak in peaks:
      start = peak - self.p.baseSamples # calculate start index
      end = peak + self.p.ccepSamples # calculate end index

      # trim indices to be within valid data range (bounds check)
      start_clipped = max(0, start)
      end_clipped = min(len(newData), end)

      data = newData[start_clipped:end_clipped] # replaced: data = newData[peak - self.p.baseSamples : peak + self.p.ccepSamples]

      # pad data to avoid size mismatch error
      if len(data) < self.p.elements:
          pad_width = self.p.elements - len(data)
          data = np.pad(data, (0, pad_width), mode='constant')

      self.computeData(data, avgPlots)

  def computeFeatures(self):
    if not hasattr (self, 'last_peak_time'):
       self.rms = np.nan
       return
    
    time_axis = self.p.x  # time axis in ms
    signal = self.data # pre-processed signal

    start_time = self.last_peak_time + 2
    end_time = 30

    mask = (time_axis >= start_time) & (time_axis <= end_time) # skip 2ms to exclude stimulation artifact
    window = signal[mask]

    if window.size == 0:
        self.rms = 0
        self.p2p = 0
    else:
        window_detrended = detrend(window, type='linear')
        self.rms = np.sqrt(np.mean(window_detrended ** 2))
        self.p2p = np.max(window_detrended) - np.min(window_detrended)

