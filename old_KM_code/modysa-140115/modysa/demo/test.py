
import io, core.configuration as C, core.geometry as G

def conformation_analysis():
    filename_pattern = '/home/murzyn/Projects/Trajectories/' + \
                       'pmb6la66-2/xtc/pmb6la66-2.part001[1-3].xtc'
    t = io.Trajectory(filename_pattern)
    print t.first, t.last
    t.limit(100000) # we start from 100ns
    # FIXME: if first in limit < t.first, t.first remains unchanged
    fxyz = t.getConfiguration()
    fcell = t.getBox()
    ftime = t.getTime()
    print ftime, fxyz.shape
    coords = C.Configuration(fxyz, fcell)
    dv = G.DistanceVector([(1,2),(2,3),(33334,33333)])
    dv(coords)
    print dv.length()
    # a molecule or a complex of molecules can to be made whole
    # according to a simple algorithm given bellow:
    #
    # 1. pass list of *n*atoms in the cluster of atoms of interest
    # 2. compute *n-1* distance_vectors between the first atom and each
    #    of the remaining atoms
    # 3. add these vectors to coordinates of the first atom
    #
    # beware: the above algorithm will fail if the largest distance within
    #         a molecule is larger than a half of the box size
    #
    # the true solution is to ensure that every valence bond in a molecule
    # is short enough; in case of a cluster of molecules, check only center of masses
    
if __name__ == '__main__': conformation_analysis()
