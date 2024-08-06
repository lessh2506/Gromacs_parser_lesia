
from PyMix import mixture
import numpy as N

def test1(FACTOR=1e5):
    fnm = '/home/murzyn/Projects/Publications/lipAH-1/dihedrals/phi__400-1199.dat'
    data = N.loadtxt(fnm)
    ds = mixture.DataSet()
    a1 = mixture.NormalDistribution(-40,10.)
    a2 = mixture.NormalDistribution(0,8.)
    m = mixture.MixtureModel(2,[.5,.5],[a1,a2])
    xx = [] # density to sample 
    for i,j in data: xx.extend([i]*int(j*FACTOR))
    ds.fromList(xx)
    m.modelInitialization(ds)
    m.EM(ds, 1000,.1)
    print m
    #clust = m.classify(ds, silent=1)
    #c1 = N.take(xx,N.nonzero(clust))
    #print N.mean(c1), N.std(c1)
    #c2 = N.take(xx,N.nonzero(clust == 0))
    #print N.mean(c2), N.std(c2)
    #return m, c1, c2

def fit_psi(FACTOR=1e5,last=-1):
    fnm = '/home/murzyn/Projects/Publications/lipAH/dihedrals/psi_400-1199-wrap.dat'
    data = N.loadtxt(fnm)[:last]
    ds = mixture.DataSet()
    a1 = mixture.NormalDistribution(172,10.)
    a2 = mixture.NormalDistribution(178,10.)
    a3 = mixture.NormalDistribution(105,10.)
    a4 = mixture.NormalDistribution(75,10.)
    a5 = mixture.NormalDistribution(75,10.)
    m = mixture.MixtureModel(5,[.43,.37,.15, .04, .01],[a1,a2,a3, a4, a5])
    xx = [] # density to sample 
    for i,j in data: xx.extend([i]*int(j*FACTOR))
    ds.fromList(xx)
    m.modelInitialization(ds)
    m.EM(ds, 1000,.1)
    print m

def fit_omega(FACTOR=1e5):
    fnm = '/home/murzyn/Projects/Publications/lipAH/dihedrals/omega_400-1199-wrap.dat'
    data = N.loadtxt(fnm)
    ds = mixture.DataSet()
    a1 = mixture.NormalDistribution(178,10.)
    a2 = mixture.NormalDistribution(68,10.)
    a3 = mixture.NormalDistribution(300,10.)
    a4 = mixture.NormalDistribution(80,10.)
    a5 = mixture.NormalDistribution(160,10.)
    a6 = mixture.NormalDistribution(30,10.)
#    a5 = mixture.NormalDistribution(70,10.)
    m = mixture.MixtureModel(6,[.4,.2,.1,.1,.1,.1],[a1,a2,a3,a4,a5,a6])
    xx = [] # density to sample 
    for i,j in data: xx.extend([i]*int(j*FACTOR))
    ds.fromList(xx)
    m.modelInitialization(ds)
    m.EM(ds, 1000,.1)
    print m
