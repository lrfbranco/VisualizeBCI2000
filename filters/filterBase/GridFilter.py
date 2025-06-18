from filters.filterBase.MasterFilter import MasterFilter
import numpy as np
import pyqtgraph as pg
from pyqtgraph.dockarea import Dock

#abstract class for any grid visualizations
class GridFilter(MasterFilter):
  def __init__(self, area, bciPath, stream):
    super().__init__(area, bciPath, stream)
  def publish(self):
    super().publish()
    self.gridPlots = pg.GraphicsLayoutWidget()

  def setConfig(self):
    super().setConfig()
    self.gridPlots.clear()
    if hasattr(self, '_maxWindows'):
      self.windows = min(self._maxWindows, self.channels)
    else:
      self.windows = self.channels
    self.windows = max(self.windows, 0)
    if self.windows == 0:
      self.numRows = 0
      self.numColumns = 0
    else:
      self.numColumns = int(np.floor(np.sqrt(self.windows)))      # this is where the grid layout is coming from (sqrt of # windows to make square)
      self.numRows = int(np.ceil(self.windows / self.numColumns))
    pass