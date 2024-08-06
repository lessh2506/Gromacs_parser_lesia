# by Marcin Kurdziel, 2009


import numpy as N
import gridwrapper

__all__=['gridsearch']

class gridsearch:
  __grid = None
  __atomnum = None

  # taken from gromacs manual
  def checkbox(self, box):
    if (box[0][1] != 0) or \
       (box[0][2] != 0) or \
       (box[1][2] != 0) or \
       (box[0][0] <= 0) or \
       (box[1][1] <= 0) or \
       (box[2][2] <= 0) or \
       (abs(box[1][0]) > 0.5*box[0][0]) or \
       (abs(box[2][0]) > 0.5*box[0][0]) or \
       (abs(box[2][1]) > 0.5*box[1][1]):
      raise ValueError, 'incorrect box'

  def buildgrid(self, box, coords, rcut):
    if (type(box).__name__ not in ['array', 'ndarray']) or \
       (type(coords).__name__ not in ['array', 'ndarray']):
      raise TypeError, 'arguments must be of numpy.ndarray or numpy.array type'
    if (N.rank(box) != 2) or (N.rank(coords) != 2):
      raise TypeError, 'arguments must be two--dimensional arrays'
    if N.shape(box) != (3, 3):
      raise ValueError, 'box must be a 3x3 array'
    (cx, cy) = N.shape(coords)
    if (cx <= 0) or (cy != 3):
      raise ValueError, 'coordinates must be a nx3 (n>0) array'
    self.checkbox(box)
    if rcut <= 0:
      raise ValueError, 'rcut must be > 0'

    if self.__grid == None:
      self.__grid = gridwrapper.py_grid()
    self.__grid.build_grid(box, coords, rcut)
    self.__atomnum = cx

  def neighbors(self, atomnum):
    if (atomnum < 0) or (atomnum >= self.__atomnum):
      raise ValueError, 'atom number is out of range'

    return self.__grid.find_neighbors(atomnum)
