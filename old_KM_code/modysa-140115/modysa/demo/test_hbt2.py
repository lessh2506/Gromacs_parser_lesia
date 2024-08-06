import modysa.io.hbonds as HBT
import modysa.io.trajectory as T
import modysa.io.index as NDX
import modysa.analyse.interactions as I
import time

def test_create(first=568000,last=568449):
    fn_xtc = '/data/T/lipAH-1a/lipAH-1a.part0040.xtc'
    demo_dir = '/home/murzyn/Projects/Software/modysa/modysa/demo/'
    fn_ndx = demo_dir + 'data/murutka-2.ndx'
    ndx = NDX.IndexFile()
    ndx.read(fn_ndx)
    hb = I.HydrogenBonds2(ndx)
    traj = T.Trajectory(fn_xtc)
    traj.limit(first,last) # [(567224.0, 576640.0)]
    block_size = 100
    hbtraj = HBT.HBTrajectory()
    hbtraj.open('%sdata/LA-1a_%d-%d.hbt' % (demo_dir,first,last),'w')

    next = True
    iblock = 0
    start_time = time.time()

    while next:
        configuration = traj.getConfiguration()
        cell = traj.getCell()
        frame_time = traj.getTime()

        hb.find(configuration,cell)
        hbtraj.add(hb.data, frame_time)
        print iblock, len(hb.data), len(hbtraj.data)
        iblock += 1

        if iblock % block_size == 0:
            print 'Speed:', (time.time()-start_time)/block_size
            start_time = time.time()
            hbtraj.sync()
            print 'Writing:', frame_time, iblock, time.time()-start_time
            start_time = time.time()
        
        next = traj.readFrame()
    if hbtraj.data: hbtraj.sync()
    hbtraj.close()

    
