from modysa.analyse import interactions as I
from modysa.io import trajectory as T
from modysa.io.gromacs import xtcio as X
from modysa.io import index as NDX
#from modysa.io import hbonds as HBT
import numpy as N, time

def test():

    #filename = "/scratch/T/lipAH-1/lipAH_000.xtc"
    filename = '/home/murzyn/Projects/Software/mdutils/lipAH_000.xtc' 
    traj = T.Trajectory(filename)
    #ndx = NDX.IndexFile('data/hb_wb-lipAH-1.ndx')
    ndx = NDX.IndexFile()
    ndx.read('data/murutka.ndx')
    hb = I.HydrogenBonds2(ndx)
    nframes = 10
    start_time = time.time()
    results = []
    for i in range(nframes):
        configuration = traj.getConfiguration()
        cell = traj.getCell()
        hb.find(configuration, cell)
        #print i, traj.getTime(), len(hb.data)
        results.append((traj.getTime(), len(hb.data)))
        traj.readFrame()
    total_time = time.time() - start_time
    print ' Done in %.2f sec.' % total_time
    print ' %.4f sec. per configuration' % (total_time/nframes,)
    a = N.array(results)
    print min(a[:,0]), max(a[:,0])
    print min(a[:,1]), max(a[:,1])
    print sum(a[:,1])/len(a)
    return hb

def test2():
    fn_xtc = '/datac1/murzyn/lipAH-1a/xtc/lipAH-1a.part0042.xtc'
    fn_ndx = 'data/murutka-2.ndx'
    ndx = NDX.IndexFile()
    ndx.read(fn_ndx)
    hb = I.HydrogenBonds2(ndx)
    traj = T.Trajectory(fn_xtc)
    nframes = 10
    start_time = time.time()
    results = []
    for i in range(nframes):
        configuration = traj.getConfiguration()
        cell = traj.getCell()
        hb.find(configuration, cell)
        #print i, traj.getTime(), len(hb.data)
        results.append((traj.getTime(), len(hb.data)))
        traj.readFrame()
    total_time = time.time() - start_time
    print ' Done in %.2f sec.' % total_time
    print ' %.4f sec. per configuration' % (total_time/nframes,)
    a = N.array(results)
    print min(a[:,0]), max(a[:,0])
    print min(a[:,1]), max(a[:,1])
    print sum(a[:,1])/len(a)



if __name__ == '__main__': 
    hb = test2()
    #print len(hb.data), ' HBonds found'
    #print hb.data
