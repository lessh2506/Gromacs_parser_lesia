
import numpy as N
import modysa.io.hbonds as HBT

def test_create():
    hbt = HBT.HBTrajectory()
    hbonds = [ [(1,2), (2,3), (3,4), (15,16)],
               [(1,2), (3,4), (18,19)],
               [(1,2), (2,3), (18,19)]]
    for i in range(len(hbonds)):
        hbt.add(hbonds[i], i)
    hbt.open('test.nc','w')
    hbt.write()
    hbt.close()
    print 'test_create PASSED'

def test_append():
    hbt = HBT.HBTrajectory()
    hbt.data = {(1,2): [[20,5],[30,18]], (32,91): [[20, 7]], 
                (33,92): [[20,8],[50,4]], (14,15): [[50,2]]}
    hbt.open('test.nc','a')
    hbt.write()
    hbt.nc.last = 54 # max for (33,92)
    hbt.close()
    print 'test_append PASSED'

def test_read():
    hbt = HBT.HBTrajectory()
    hbt.open('test.nc','r')
    print hbt.pairs
    print hbt.nc.variables['pairs'][:]
    print hbt.nc.variables['events'][:]
    hbt.close()
    print 'test_read PASSED'

def test_misc():
    hbt = HBT.HBTrajectory()
    hbt.open('test.nc','r')
    print hbt.trajectory(1,2)
    print hbt.events(33,92)
    pairs = ((1,2),(33,92))
    print hbt.number(pairs,first=10)
    print hbt.number(pairs)
    print hbt.number(pairs,last=45)
    print hbt.number(pairs,step=2)
    hbt.close()
    print 'test_misc PASSED'

def clean():
    import os
    os.remove('test.nc')

if __name__ == '__main__': 
    clean()
    test_create()
    test_append()
    test_read()
    test_misc()

