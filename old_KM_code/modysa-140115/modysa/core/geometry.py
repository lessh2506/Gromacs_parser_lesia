
import numpy as N

class Vector:

    def __init__(self):
        self.array = None

    def __call__(self, xyz):
        self.array = N.array(xyz)
        assert self.array.shape[-1] == 3

    def length(self):
        assert self.array is not None
        return N.sqrt(N.add.reduce(self.array*self.array,-1))

    def projection(self, vector):
        assert isinstance(vector, Vector) and self.array is not None
        shape = self.array.shape
        array = self.array * vector.array[N.newaxis,:]
        self.array = array.reshape(shape)

    def cross(self, other):
        assert isinstance(other, Vector)
        v = Vector()
        v(N.cross(self.array, other.array))
        return v

    def angle2(self, other):
        assert isinstance(other, Vector)
        v1 = self.array
        v2 = other.array
        assert v1.shape[-1] == v2.shape[-1]
        c = N.cross(v1,v2)
        # dot perp prod in 2d calculates z-component of the normal vector
        cross_magn = N.sqrt(N.add.reduce(c*c,-1)) 
        dot_prod = N.add.reduce(v1*v2,-1)
        angle = N.arctan2(cross_magn,dot_prod)
        # which sign?
        # it seems we have to have reference vector.. right handed system
        # see Answer 1: http://www.techques.com/question/1-1180849/Signed-Angle-in-3D-Vectors
        reference = other.rotate(self,N.pi/2)
        sign = N.add.reduce(reference*c, -1) # projection
        # FIXME: self.projection()
        angle = N.where(sign < 0, -angle, angle)
        return angle

    def normalize(self):
        length = self.length()
        assert N.all(length != 0)
        self.array = self.array/length[...,N.newaxis]
        #if self.array.ndim == 1:
        #    self.array = self.array/length
        
            

    def rotate(self, axis, theta):
        # theta in Rad
        # http://stackoverflow.com/a/6802723
        # http://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_parameters
        assert isinstance(axis, Vector)
        axis.normalize()
        a = N.cos(theta/2.)
        bcd = -axis.array*N.sin(theta/2.)
        b, c, d = map(N.array,zip(*bcd[N.newaxis,...]))
        rotmatrix = N.array([[a*a+b*b-c*c-d*d, 2.*(b*c-a*d), 2*(b*d+a*c)],
                             [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                             [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
        rotvec = Vector()
        #rotvec(N.dot(rotmatrix, self.array))
        # FIXME: problem with multiplication, shapes not aligned
        # N.dot(rotmatrix[...,0],nv2.array[0])
        self.rotmatrix = rotmatrix
        return rotvec

    def angle(self, other, reference=(0,0,1)):
        assert isinstance(other, Vector)
        v1 = self.array
        v2 = other.array
        if isinstance(reference, Vector): vn = N.array(reference.array)
        else: vn = N.array(reference)
        if len(vn.shape) == 1: 
            vn = vn[N.newaxis,:] # one reference vector only
        assert v1.shape[-1] == v2.shape[-1]
        if v1.shape[-1] > 2: signed = True
        else: signed = False
        denom = N.add.reduce(v1*v1,-1) * N.add.reduce(v2*v2,-1)
        # singularity? A warning from numpy:
        # Warning: invalid value encountered in divide
        denom = N.where(denom == 0, 1, denom)
        cosa = N.add.reduce(v1*v2,-1) / N.sqrt(denom) 
        cosa = N.where(N.less(N.where(N.greater(cosa,-1.),cosa,-1.),1.),
                       cosa,1.)
        angle = N.arccos(cosa)
        # http://stackoverflow.com/a/5190354
        if signed:
            c = self.cross(other).array
            sign = N.add.reduce(vn*c, -1) # dot product
            # if N.any(sign == 0): 
            #     err = N.nonzero(N.where(sign == 0, 1,0))
            #     if N.any(c[err]) != 0:
            #         # projection
            #         print 'Warning: possibly wrong reference vector'
            #         print c[err], v1[err], v2[err], angle[err]
            #     # else: anti-parallel, i.e. pi/-pi which doesn't matter
            angle = N.where(sign < 0, -angle, angle)
        return angle
    

class DistanceVector(Vector):

    def __init__(self, pairs):
        self.setA, self.setB = zip(*pairs)
        self.array = None

    def __len__(self):
        return self.array.shape[0]

    def __call__(self, configuration):
        """ configuration is an instance of 
            modysa.core.configuration.Configuration object
        """
        xyz1 = N.take(configuration.array, self.setA, axis=0)
        xyz2 = N.take(configuration.array, self.setB, axis=0)
        self.array = xyz2 - xyz1
        # handling PBC
        if configuration.cell is not None:
            cell = configuration.cell.parameters()
            invert = configuration.cell.invert()
            df = N.zeros(self.array.shape, N.float)
            for i in range(3):
                df[...,i] = N.add.reduce(invert[i*3:(i+1)*3] * self.array, -1)
                unwrapped = True
                while unwrapped:
                    df[...,i] = N.where(N.greater(df[...,i],.5),
                                        df[...,i]-1., df[...,i])
                    df[...,i] = N.where(N.less_equal(df[...,i],-.5),
                                        df[...,i]+1.,df[...,i])
                    unwrapped = N.add.reduce(N.ravel(\
                            N.logical_and(N.greater(df[...,i],.5),
                                          N.less_equal(df[...,i],-.5))))
            for i in range(3): 
                self.array[...,i] = N.add.reduce(cell[i*3:(i+1)*3]*df,-1)

    def length(self):
        assert self.array is not None
        return N.sqrt(N.sum(self.array**2,-1))


class Rotation:

    def __init__(self):
        self.matrix = None

    def fromQuaternions(self, quaternions):
        """ make the rotation matrix from quaternions, see
        core.stack.RigidBodyTrajectory.rotationMatrix """
        # if len(quaterions.shape) == 2: set rotation *matrices* for each
        # quaternion provided, the same applies to other from... methods
        self.matrix = None

    def fromRotationMatrix(self, matrix):
        self.matrix = matrix

    def fromAxisAndAngle(self, axes=None, thetas=None):
        # see below (setMatrix)
        pass

    def setMatrix(self, axes=None, thetas=None, array=None):
        if axes is not None and thetas is not None:
            assert isinstance(axes, Vector)
            axes.normalize()
            a = N.cos(thetas/2.)
            bcd = -axes.array*N.sin(thetas/2.)
            b, c, d = map(N.array,zip(*bcd[N.newaxis,...]))
            matrix = N.array([[a*a+b*b-c*c-d*d, 2.*(b*c-a*d), 2*(b*d+a*c)],
                              [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                              [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
            self.matrix = N.swapaxes(matrix, 0, 2)
            #self.matrix = N.swapaxes(matrix, 1, 2)
        else:
            assert array is not None
            self.matrix = array

    def __call__(self, vectors):
        """ apply rotation defined by rotation matrix (self.matrix)
        to provided 3d vectors """
        # FIXME: it should be possible to rotate n vectors with n matrices
        # or rotate n vectors with just one matrix, n might be 1,
        # assert self.matrix.shape[0] == vectors.shape[0]
        assert self.matrix is not None
        # check this with core.stack.Stack.superpose: 
        #    dot(inverse(D),transpose(vectors))
        return N.add.reduce(self.matrix*vectors.array, -1) # dot product

def angle(vector1,vector2):
    """ angle [rad] between two vectors """

    denom = N.add.reduce(vector1*vector1,-1) * N.add.reduce(vector2*vector2,-1)
    # singularity? A warning from numpy:
    # Warning: invalid value encountered in divide
    denom = N.where(denom == 0, 1, denom)
    cosa = N.add.reduce(vector1*vector2,-1) / N.sqrt(denom) 
    cosa = N.where(N.less(N.where(N.greater(cosa,-1.),cosa,-1.),1.),cosa,1.)

    return N.arccos(cosa)

def dihedral( vector1, vector2, vector3):
    
    #a  = cross(vector1, vector2)
    #b  = cross(vector2, vector3)
    a = N.cross(vector1, vector2)
    b = N.cross(vector2, vector3)
    rarb = N.sqrt(N.add.reduce(a*a,-1)*N.add.reduce(b*b,-1))
    cos  = N.add.reduce(a*b,-1)/rarb
    cos  = N.where(N.greater(-1.,cos),-1.,cos)
    cos  = N.where(N.less(1.,cos),1.,cos)
    angle = N.arccos(cos)
    sign  = N.add.reduce(vector1*b,-1)
    angle = N.where(N.less(0.,sign),angle,-angle)
    return angle

def cross(v1,v2):
    """ cross product of two vectors """
    a = N.zeros((max(len(v1),len(v2)),3), N.float)
    a[:,0] = v1[:,1]*v2[:,2] - v1[:,2]*v2[:,1]
    a[:,1] = v1[:,2]*v2[:,0] - v1[:,0]*v2[:,2]
    a[:,2] = v1[:,0]*v2[:,1] - v1[:,1]*v2[:,0]
    return a


if __name__ == '__main__': 
    frame = N.array([[ 56.7  ,   0.25 ,  16.95 ],
                     [ 60.45 ,  -2.048,  18.75 ],
                     [ 58.08 ,   2.2  ,  17.7  ],
                     [ -17.092,  -0.81 ,  92.523  ],
                     [ 59.44 ,  -0.678,  20.54 ],
                     [ 55.57 ,   1.67 ,  20.03 ],
                     [ 57.77 ,   0.79 ,  20.86 ],
                     [ 57.02 ,  73.61 ,  17.65 ],
                     [ 58.1  ,  75.42 ,  18.6  ]])
    cell = N.array([[ 73.252,   0.   ,   0.   ],
                    [  0.   ,  76.258,   0.   ],
                    [  0.   ,   0.   ,  72.723]])
    import modysa.core.configuration as C
    configuration = C.Configuration(frame, cell)
    dv = DistanceVector(zip(range(8),[8]*8)) # eight first vectors vs last one
    dv(configuration)
    print dv.length()
    

#     import pymacs as P
#     from MMTK import Units
#     fp = P.openXTC('../pe81cTc_350.xtc')
#     frame = P.readXTCFrame(fp)
#     conf = Configuration(frame['x'], cell_shape = frame['box'])
#     setA = [34827] # oxygen
#     setB = [34828, 34829] # hydrogens
#     d = Distance(setA,setB)
#     d(conf)
#     print d.array[0][0], d.array[0][1]
#     d(conf,1)
#     v1 = d.array[0][0]
#     v2 = d.array[0][1]
#     print angle(v1,v2)/Units.deg
