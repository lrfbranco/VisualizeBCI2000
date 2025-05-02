# UDPCommunication.py

import socket
import numpy as np
from PyQt5.QtCore import QObject, pyqtSignal, QThread
from dataThreads.AbstractClasses import AbstractCommunication, AbstractWorker, AbstractDataThread

class UDP(AbstractCommunication):
    def __init__(self, host, port):
        # instantiate the data thread and worker
        self._acqThr = UDPDataThread()
        self._worker = UDPWorker(host, port)

    @property
    def worker(self):
        return self._worker

    @property
    def acqThr(self):
        return self._acqThr

    def evaluate(self, string):
        # no real “state” to query over UDP; stubbed
        return False

class UDPWorker(AbstractWorker):
    def __init__(self, host, port):
        super().__init__()
        self.host = host
        self.port = port
        self._isRunning = True

    def run(self):
        # immediately tell the data thread where to listen
        addr = f"{self.host}:{self.port}"
        self.initSignal.emit(addr)
        # idle until stopped
        while self._isRunning:
            QThread.msleep(100)
        self.disconnected.emit()

    def stop(self):
        self._isRunning = False

class UDPDataThread(AbstractDataThread):
    def __init__(self):
        super().__init__()
        self._isRunning = True
        self._props_set = False

    def initialize(self, address):
        host, port = address
        self.s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.s.bind((host, int(port)))
        self.s.settimeout(0.5)
        self.printSignal.emit(f"Listening for UDP on {host}:{port}")

    def run(self):
        while self._isRunning:
            try:
                packet, _ = self.s.recvfrom(65536)
                text = packet.decode('utf-8')
                # PROPS:ch1,ch2,ch3
                if text.startswith('PROPS:'):
                    names = text[len('PROPS:'):].split(',')
                    # emit 1 “element” (sample), with these channel names
                    self.propertiesSignal.emit(1, names)
                    self._props_set = True

                # DATA:v1,v2,v3
                elif text.startswith('DATA:') and self._props_set:
                    vals = [float(x) for x in text[len('DATA:'):].split(',')]
                    arr = np.array(vals, dtype=float).reshape((len(vals), 1))
                    self.dataSignal.emit(arr)

                else:
                    self.printSignal.emit(f"Unknown message: {text}")

            except socket.timeout:
                continue
            except Exception as e:
                self.printSignal.emit(f"Error: {e}")

    def stop(self):
        self._isRunning = False
        self.s.close()
