import ctypes as C
import numpy as N
import gmxiobase as G

__all__=['TRNReader']

class TRNHeader(C.Structure):

  _fields_ = [("ir_size", C.c_int), \
              ("e_size", C.c_int), \
              ("box_size", C.c_int), \
              ("vir_size", C.c_int), \
              ("pres_size", C.c_int), \
              ("top_size", C.c_int), \
              ("sym_size", C.c_int), \
              ("x_size", C.c_int), \
              ("v_size", C.c_int), \
              ("f_size", C.c_int), \
              ("natoms", C.c_int), \
              ("step", C.c_int), \
              ("nre", C.c_int), \
              ("t", C.c_float), \
              ("l", C.c_float)]

class TRNReader:

  def __init__(self):
    self.__fileDescr = -1

    self.__header = TRNHeader()
    self.__box = G.Matrix()
    #self.__vectorArrayT = C.c_float
    self.__coords = G.Vector()
    self.__veloc = G.Vector()
    self.__forces = G.Vector()
    self.__frameStatus = G.Boolean(-1)

    self.__header.natoms = -1


  def openTRN(self, filename):

    if not os.path.exists(filename):
      raise IOError, 'Could not open TRJ/TRR file - '+filename

    if self.__fileDescr != -1: self.closeTRN()

    self.__fileDescr = Libgmx.open_trn(filename, "r")
    if self.__fileDescr < 0:
      raise IOError, 'Could not open TRJ/TRR file - '+filename
      self.__fileDescr = -1

    self.readFrame()

  def closeTRN(self):

    if self.__fileDescr == -1:
      raise IOError, 'No TRJ/TRR file opened'

    Libgmx.close_trn(self.__fileDescr)
    self.__fileDescr = -1

    self.__header.natoms = -1
    self.__header.step = -1
    self.__header.t = -1
    self.__header.l = -1
    self.__frameStatus.value = -1


  def readFrame(self):

    if self.__fileDescr == -1:
      raise IOError, 'No TRJ/TRR file opened'

    res = Libgmx.fread_trnheader(self.__fileDescr, 
                                 C.pointer(self.__header), 
                                 C.pointer(self.__frameStatus))
    if res == 0:
      raise IOError, 'Error while reading TRJ/TRR frame header'

    if self.__frameStatus.value == 0:
      raise IOError, 'Incorrect TRJ/TRR frame header'

    self.__vectorArrayT = (c_float * self.__header.natoms)*3
    self.__coords = self.__vectorArrayT()    
    self.__veloc = self.__vectorArrayT()    
    self.__forces = self.__vectorArrayT()    

    res = Libgmx.fread_htrn(self.__fileDescr, 
                            C.pointer(self.__header), C/pointer(self.__box),
                            C.pointer(self.__coords), C.pointer(self.__veloc), 
                            C.pointer(self.__forces))

    if res == 0:
      raise IOError, 'Error while reading TRJ/TRR frame'

  def getAtomsNum(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return self.__header.natoms


  def getStep(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return self.__header.step


  def getTime(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return self.__header.t


  def getLambda(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return self.__header.l


  def getFrameStatus(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return self.__frameStatus.value


  def getBox(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return array(self.__box, Float)


  def getCoords(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return reshape(array(self.__coords, Float), (self.__header.natoms, 3))


  def getVeloc(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return reshape(array(self.__veloc, Float), (self.__header.natoms, 3))


  def getForces(self):
    if self.__header.natoms == -1:
      readFrame(self)
    return reshape(array(self.__forces, Float), (self.__header.natoms, 3))
