import numpy as N
import modysa.misc as M

U = M.Units

class Cell:

    __defaults = {'a': None, 'b': None, 'c': None,
                  'bc': 90., 'ac': 90., 'ab': 90.}

    def __init__(self):
        self.vectors = N.zeros((3,3),N.float)

    def setFromBasisVectors(self, vectors):
        assert type(vectors).__name__ == 'ndarray' and vectors.shape == (3,3)
        self.vectors = vectors

    def setFromLatticeConstants(self, **params):
        for k in self.__defaults.keys():
            setattr(self, k, M.setDefault(params,k,self.__defaults[k]))
        if self.a is None or self.b is None or self.c is None:
            return
        # see nMoldyn3.0 Mathematics.Geometry.basisVectors
        self.vectors[0] = [self.a,0,0]
        self.vectors[1] = [N.cos(self.ab*U.deg)*self.b, 
                           N.sin(self.ab*U.deg)*self.b, 0]
        cx = N.cos(self.ac*U.deg)
        cy = (N.cos(self.bc*U.deg) - \
                  N.cos(self.ac*U.deg)*N.cos(self.ab*U.deg)) / \
                  N.sin(self.ab*U.deg)
        cz = N.sqrt(1.0 - cx**2 - cy**2)
        self.vectors[2] = self.c * N.array([cx,cy,cz])

    def volume(self):
        assert N.all(N.any(self.vectors,axis=0)) # no 0,0,0 vectors
        x,y,z = N.vectors
        # triple product
        return N.dot(N.cross(x,y), z)

    def vectorArea(self, norm=None):
        #  temporary solution (works on selection of two basis vectors
        #  out of three)
        # FIXME; 
        assert N.all(N.any(self.vectors,axis=0)) # no 0,0,0 vectors
        # 1. which two vectors go to cross product
        a, b = N.take(self.vectors, N.nonzero(N.logical_not(norm))[0], axis=0)
        c = N.cross(a,b)
        # 2. area is just the length of cross product vector
        return N.sqrt(N.sum(c*c))

    def parameters(self):
        return N.ravel(N.transpose(self.vectors))
        #return N.ravel(N.transpose([list(s) for s in self.vectors]))
        
    def invert(self):
        """
        Calculate determinant and inverse of the coordinate transformation
        matrix for parallelepipedic universes.
        """
        cell = self.parameters()
        assert cell is not None

        data = N.zeros((9,), N.float)
        data[0] = cell[4]*cell[8] - cell[7]*cell[5]
        data[3] = cell[6]*cell[5] - cell[3]*cell[8]
        data[6] = cell[3]*cell[7] - cell[6]*cell[4]
        data[1] = cell[7]*cell[2] - cell[1]*cell[8]
        data[4] = cell[0]*cell[8] - cell[6]*cell[2]
        data[7] = cell[6]*cell[1] - cell[0]*cell[7]
        data[2] = cell[1]*cell[5] - cell[4]*cell[2]
        data[5] = cell[3]*cell[2] - cell[0]*cell[5]
        data[8] = cell[0]*cell[4] - cell[3]*cell[1]
        #
        det = cell[0]*data[0] + cell[1]*data[3] + cell[2]*data[6]
        if abs(det) > 0.: r = 1./det
        else: r = 0.
        data = data * r
        return data


class Configuration:

    def __init__(self, xyz, cell = None):
        assert isinstance(cell, Cell)
        self.update(xyz, cell)

    def update(self, xyz, cell=None):
        self.array = N.array(xyz)
        self.cell = cell 
        # FIXME: the previously written code 
        #    needs reverting to the present version of Configuration/Cell
        #      self.cell.vectors <- self.cell_vectors
        #      self.cell.parameters() = self.cell_parameters

