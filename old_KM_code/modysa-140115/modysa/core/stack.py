import numpy as N
import LinearAlgebra
# import Scientific.Geometry as G # Vector, Tensor 


version = '2008.03.19'

class RigidBodyTrajectory:

    def __init__(self,rms,cms,quaternions,reference):

        self.rms = rms
        self.cms = cms
        self.quaternions = quaternions
        self.reference = reference

    def rotationMatrix(self):

        q = self.quaternions
        D = N.zeros((len(q),3,3),N.Float)
    
        D[:,0,0] = -2.*q[:,2]**2 - 2.*q[:,3]**2 + 1.
        D[:,0,1] = 2.*(-q[:,0]*q[:,3] + q[:,1]*q[:,2])
        D[:,0,2] = 2.*(q[:,0]*q[:,2] + q[:,1]*q[:,3])

        D[:,1,0] = 2.*(q[:,0]*q[:,3] + q[:,1]*q[:,2])
        D[:,1,1] = -2.*q[:,1]**2 - 2.*q[:,3]**2 + 1.
        D[:,1,2] = 2.*(-q[:,0]*q[:,1] + q[:,2]*q[:,3])

        D[:,2,0] = 2.*(-q[:,0]*q[:,2] + q[:,1]*q[:,3])
        D[:,2,1] = 2.*(q[:,0]*q[:,1] + q[:,2]*q[:,3])
        D[:,2,2] = -2.*q[:,1]**2 - 2.*q[:,2]**2 + 1.

        return D

class Stack:

    def __init__(self, trajectory, mass=None):
        """ atom 
                contains xyz coordinates for nat atoms 
                from nmol molecules
            mass
                contains atomic masses for corresponding atoms
                if None, all masses are set to 1 
        """
        self.trajectory = trajectory

        if mass is None: 
            self.mass = N.ones((self.trajectory.shape[1],),
                               N.Float)
        else: 
            assert(len(mass) == self.trajectory.shape[1])
            self.mass = N.array(mass,N.Float)

        self.rbt = None

    def getRigidBodyTrajectory(self, reference=None):
        """ reference 
                is an array of xyz coordinates for nat atoms 
                of a reference molecule/atom cluster 
        """
        if reference is None: reference = self.trajectory[0]

        weights = self.mass/N.add.reduce(self.mass)
        ref_cms = N.add.reduce(weights[:,N.NewAxis]*reference)
        traj_cms = N.add.reduce(weights[N.NewAxis,:,N.NewAxis]*\
                                    self.trajectory,1)
        r_ref = reference - ref_cms
        r_traj = self.trajectory - traj_cms[:,N.NewAxis,:]

        steps = self.trajectory.shape[0]
        possq = N.add.reduce(weights*N.add.reduce(r_traj*r_traj,-1),-1)+\
                N.add.reduce(weights*N.add.reduce(r_ref[N.NewAxis,:]*\
                                        r_ref[N.NewAxis,:],-1),-1)
        cross = N.add.reduce(weights[N.NewAxis,:,N.NewAxis,N.NewAxis]*\
                                 r_traj[:,:,:,N.NewAxis]*\
                                 r_ref[N.NewAxis,:,N.NewAxis,:],1)

        k = N.zeros((steps, 4, 4), N.Float)
        k[:, 0, 0] = -cross[:, 0, 0]-cross[:, 1, 1]-cross[:, 2, 2]
        k[:, 0, 1] = cross[:, 1, 2]-cross[:, 2, 1]
        k[:, 0, 2] = cross[:, 2, 0]-cross[:, 0, 2]
        k[:, 0, 3] = cross[:, 0, 1]-cross[:, 1, 0]
        k[:, 1, 1] = -cross[:, 0, 0]+cross[:, 1, 1]+cross[:, 2, 2]
        k[:, 1, 2] = -cross[:, 0, 1]-cross[:, 1, 0]
        k[:, 1, 3] = -cross[:, 0, 2]-cross[:, 2, 0]
        k[:, 2, 2] = cross[:, 0, 0]-cross[:, 1, 1]+cross[:, 2, 2]
        k[:, 2, 3] = -cross[:, 1, 2]-cross[:, 2, 1]
        k[:, 3, 3] = cross[:, 0, 0]+cross[:, 1, 1]-cross[:, 2, 2]

        for i in range(1, 4):
            for j in range(i): k[:, i, j] = k[:, j, i]
        N.multiply(k, 2., k)
        for i in range(4): N.add(k[:,i,i], possq, k[:,i,i])

        quaternions = N.zeros((steps, 4), N.Float)
        fit = N.zeros((steps,), N.Float)
        for i in range(steps):
            e, v = LinearAlgebra.eigenvectors(k[i])
            j = N.argmin(e)
            if e[j] < 0.: fit[i] = 0.
            else: fit[i] = N.sqrt(e[j])
            if v[j,0] < 0.: quaternions[i] = -v[j]
            else: quaternions[i] = v[j]

        self.rbt = RigidBodyTrajectory(fit,traj_cms,quaternions,reference)

    def superpose(self, imol, conf=None):
        """ superpose imolth configuration on the reference structure or
            if conf is given, imolth rotation/translation transformation
            is applied to it
        """
        if self.rbt is None: return
        if conf is None: xyz = self.trajectory[imol]
        else: xyz = conf # arrays
        D = self.rbt.rotationMatrix()
        gj = N.dot(LinearAlgebra.inverse(D[imol]),
                   N.transpose(xyz-self.rbt.cms[imol]))
        weights = self.mass/N.add.reduce(self.mass)
        # add cms of reference
        xyz = N.transpose(gj) + \
            N.add.reduce(weights[:,N.NewAxis]*self.rbt.reference)
        return xyz

    def rmsdPerAtom(self, imol, reference=None):
        """  returns global RMSD and RMSD per Atom
        """
        if self.rbt is None: return
        if reference is None: reference = self.rbt.reference
        r = self.trajectory[imol] - reference
        # weights: max is 1.0
        dr = N.add.reduce(r*r,-1)*self.mass/N.add.reduce(self.mass)*\
            len(self.mass)
        rmsd = N.sqrt(N.add.reduce(dr)/len(dr))
        return rmsd, dr

class RadialDistributionFunction:

    def __init__(self,confA,confB=None,r_max=None,width=0.02,cell=None):
        """ cell is Tensor - three vectors describing unit cell
        """
        # d = distance(configuration_array, box_tensor, atom_pairs=None)
        r = distance(confA, confB) # , cell) # PBC to be implemented
        bins = N.arange(0.,r_max,width)
        # histogram
        n = N.searchsorted(N.sort(r),bins)
        n = N.concatenate([n,[len(r)]])
        rdf = (n[1:] - n[:-1])[:-1]
        # normalize RDF
        pairs = len(r)
        if cell is not None and hasattr(cell,'is_tensor'):
            vx, vy, vz = cell
            volume = vx.asVector()*vy.asVector().cross(vz.asVector())
        else: volume = 4/3.*N.pi*r_max**3 # hmm
        factor = 4./3. * N.pi * pairs * (bins[1:]**3 - bins[:-1]**3)/volume
        self.data = rdf/factor

# ---



# ---

if __name__ == '__main__':

    print ' Stack module ', version
    print ' Testing superposition...'

    from Scientific.IO import PDB
    from Scientific import Statistics

    fA, fB = '1Z1Y.pdb', '1Z27.pdb'
    pdb_1Z1Y = PDB.Structure(fA)
    pdb_1Z27 = PDB.Structure(fB)

    pA = pdb_1Z1Y.residues[20:30]
    pB = pdb_1Z27.residues[20:30]

    tmol_A = N.array([r.atoms['CA'].position.array for r in pA])
    tmol_B = N.array([r.atoms['CA'].position.array for r in pB])

    print ' . Sequence %s' % fA,".".join([r.name for r in pA])
    print ' . Sequence %s' % fB,".".join([r.name for r in pB])

    trajectory = N.array([tmol_A,tmol_B])
    print ' Number of molecules in the Stack:',len(trajectory)
    print ' Number of atoms in each of molecules:', trajectory.shape[1]
    
    s = Stack(trajectory,mass=N.array([12.]*10))
    s.getRigidBodyTrajectory()
    print ' Center of mass before superposition'
    print ' . %s (reference):' % fA, s.rbt.cms[0]
    print ' . %s (fitted):'    % fB, s.rbt.cms[1]
    print ' RMSD after superposition (from quaternions): ', s.rbt.rms[-1]

    rmsd, rmsd_at = s.rmsdPerAtom(0,s.trajectory[1])
    print ' The total (average) RMSD per atom before superposition: %f (%f)' %\
        (rmsd, Statistics.mean(rmsd_at))
    rmsd, rmsd_at = s.rmsdPerAtom(0,s.superpose(1))
    print ' The total (average) RMSD after superposition: %f (%f)' %\
        (rmsd, Statistics.mean(rmsd_at))
