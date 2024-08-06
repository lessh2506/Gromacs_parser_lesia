
import glob

def getFilenames(pattern='/data1/Trajectories/murzyn/lipAH-2/run/xtc/*.xtc'):
    fnms = glob.glob(pattern)
    fnms.sort()
    return fnms
