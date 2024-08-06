import ctypes as C
import gmxiobase as G
import numpy as N
import os.path, re
import modysa.misc as M

__all__ = ['XTCReader']

class XTCReader:

  def __init__(self):

    self.__fileDescr = -1

    self.__atomsNum = C.c_int(-1)
    self.__step = C.c_int(-1)
    self.__time = C.c_float(-1)
    self.__box = G.Matrix()
    self.__coordsArrayT = C.c_float
    self.__coords = G.PVector()
    self.__prec = C.c_float(-1)
    self.__frameStatus = G.Boolean(-1)


  def open(self, filename):

    if not os.path.exists(filename):
      raise IOError, 'Could not open XTC file:' + filename
    if self.__fileDescr != -1: self.close()
    self.__fileDescr = G.Libgmx.open_xtc(filename, "r")
    if self.__fileDescr < 0:
      self.__fileDescr = -1
      raise IOError, 'Could not open XTC file:' + filename
    self.readFrame()

  def close(self):

    if self.__fileDescr == -1:
      raise IOError, 'No XTC file currently opened'
    G.Libgmx.close_xtc(self.__fileDescr)
    self.__fileDescr = -1

    if self.__atomsNum.value != -1: G.LibC.free(self.__coords)

    self.__atomsNum.value = -1
    self.__step.value = -1
    self.__time.value = -1
    self.__prec.value = -1
    self.__frameStatus.value = -1


  def readFrame(self):

    if self.__fileDescr == -1:
      raise IOError, 'No XTC file opened'

    if self.__atomsNum.value == -1:
      if G.Libgmx.read_first_xtc(self.__fileDescr, C.pointer(self.__atomsNum), 
                                 C.pointer(self.__step), 
                                 C.pointer(self.__time), self.__box, 
                                 C.pointer(self.__coords), 
                                 C.pointer(self.__prec), 
                                 C.pointer(self.__frameStatus)) == 0:
        raise IOError, 'Error while reading XTC frame'
      self.__coordsArrayT = (C.c_float * self.__atomsNum.value) * 3
    else:
      if G.Libgmx.read_next_xtc(self.__fileDescr, self.__atomsNum, 
                                C.pointer(self.__step), 
                                C.pointer(self.__time),
                                self.__box, self.__coords, 
                                C.pointer(self.__prec), 
                                C.pointer(self.__frameStatus)) == 0:
        raise IOError, 'Error while reading XTC frame'

  def tell(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    # XTC file structure is described in libxdrf.c 
    # (just before xtc_at_header_start function);
    # be aware that consecutive frames may have different size in bytes;
    # to allow moving along the XTC trajectory one needs to record the
    # begining position of each frame in external file (e.g. traj extension)
    # and use these positions/offsets with Libgmx.gmx_fio_seek(fio, offset)
    fio = self.__fileDescr
    G.Libgmx.gmx_fio_ftell.restype = c_long
    return G.Libgmx.gmx_fio_ftell(fio)

  #   def seek(self, offset):
  #     # beware that if the header of a frame doesn't match MAGIC number
  #     # (0x000007CB or 1995) Libgmx function quits with an error msg
  #     fio = self.__fileDescr
  #     # offset is off_t which type is defined by configure (does your system support large files or not)
  #     Libgmx.gmx_fio_seek(fio, c_long(offset))

  def seekByStep(self, step):
    if self.__fileDescr < 0: return # no trajectory file opened!
    # step of MD, beware of the fact that trajectory is dumped every nth step
    fio = self.__fileDescr
    natoms = self.getAtomsNum()
    G.Libgmx.xtc_seek_frame(step, fio, natoms)
    # FIXME: what if step < 0 or step > max(traj.steps)?
    self.readFrame()

  def seekByTime(self, t):
    if self.__fileDescr < 0: return # no trajectory file opened!
    fio = self.__fileDescr
    natoms = self.getAtomsNum()
    G.Libgmx.xtc_seek_time(C.c_float(t), fio, natoms)
    # FIXME: what if t < 0 or t > max(traj.time)?
    self.readFrame()

  def getLastFrameTime(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    fio = self.__fileDescr
    natoms = self.getAtomsNum()
    # the trick is to *cast* pointers so we don't need to 
    # define FILE/xdrs rather intricate structures and 
    # it's enough to pass void pointer
    G.Libgmx.gmx_fio_getfp.restype = C.POINTER(C.c_void_p)
    FILE = G.Libgmx.gmx_fio_getfp(fio)
    #
    G.Libgmx.gmx_fio_getxdr.restype = C.POINTER(C.c_void_p)
    xdrs = G.Libgmx.gmx_fio_getxdr(fio)
    #
    if re.match('4\.[0-4]\.[0-9]+', G.gmx_version):
      func = G.Libgmx.xtc_get_last_frame_time
    elif re.match('4\.5\.[0-9]+', G.gmx_version):
      func = G.Libgmx.xdr_xtc_get_last_frame_time
    else: # guessing...
      func = G.Libgmx.xdr_xtc_get_last_frame_time
    #
    func.restype = C.c_float
    bOK = G.Boolean(0) # FIXME: check what bOK really is
    t = func(FILE, xdrs, natoms, C.pointer(bOK))
    return M.floatSignificantNumbers(t)


  def getAtomsNum(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()
    return self.__atomsNum.value


  def getStep(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()
    return self.__step.value


  def getTime(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()
    return M.floatSignificantNumbers(self.__time.value)


  def getPrec(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()
    return self.__prec.value


  def getFrameStatus(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()
    return self.__frameStatus.value


  def getCell(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()
    return N.array(self.__box, N.float)


  def getCoords(self):
    if self.__fileDescr < 0: return # no trajectory file opened!
    if self.__atomsNum.value == -1: self.readFrame()

    coordsArray = self.__coordsArrayT()
    G.LibC.memcpy(coordsArray, self.__coords, C.sizeof(self.__coordsArrayT))
    # http://stackoverflow.com/questions/4964101/pep-3118-warning-when-using-ctypes-array-as-numpy-array
    # xtcio.py:179: RuntimeWarning: Item size computed from the PEP 3118 buffer format string does not match the actual item size.
    return N.reshape(N.array(coordsArray, N.float), (self.__atomsNum.value, 3))


