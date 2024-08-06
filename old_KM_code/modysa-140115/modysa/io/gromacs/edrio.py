
import ctypes as C
import gmxiobase as G
import numpy as N
import os.path

filename = '/home/murzyn/Projects/Students/MalgorzataBerg/Ene/procarp-3.part0076.edr'

class edrFrame(C.Structure):

  _fields_ = [("t", C.c_float),
              ("step", C.c_ulong),
              ("nsteps", C.c_ulong),
              ("dt", C.c_float),
              ("nsum", C.c_int),
              ("nre", C.c_int),
              ("e_size", C.c_int),
              ("e_alloc", C.c_int),
              ("ener", C.c_void_p),
              ("nblock", C.c_int),
              ("block", C.c_void_p),
              ("nblock_alloc", C.c_int)]

## missing: enxio.h:t_enxblock, t_enxsubblock types/energy.h:t_energy

class EDRReader:

    def __init__(self):
        self.__fileDescr = -1
        self.__data = edrFrame()

    def open(self, filename):
        if not os.path.exists(filename):
            raise IOError, 'Could not open EDR file:', filename
        self.__fileDescr = G.Libgmx.open_enx(filename, "r")
        if self.__fileDescr < 0:
            self.__fileDescr = -1
            raise IOError, 'Could not open EDR file:', filename

    def close(self):
        if self.__fileDescr == -1:
            raise IOError, 'No EDR file currently opened'
        G.Libgmx.close_xtc(self.__fileDescr)
        self.__fileDescr = -1

    def read(self):
        G.Libgmx.do_enx(self.__fileDescr, C.pointer(self.__data))

x = EDRReader()
x.open(filename)
x.read()
x.close()

