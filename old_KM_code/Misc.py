#from Scientific.IO import PDB as P


def uniqueCombinations(items, n):
    """ Generators for calculating the combinations and selections
    of a number of elements from a sequence

    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465 """
    
    if n==0: yield []
    else:
        for i in xrange(len(items)-n+1):
            for cc in uniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc


def lj(distance,sigma,epsilon):
    sr = (sigma/distance)**6
    return 4*epsilon*(sr**2 - sr)

def getCoulombEnergy(distance,qi,qj):
    return 138.935485*qi*qj/distance

def origE(configuration,qch2=0.005,qch3=-.07):
    pass

def fitE(configuration):
    pass



class PotentialEnergy:

    def __init__(self,pdbfilename,system,topology):
        pdbfile = P.Structure(pdbfilename)
        self.system = []
#        for residue in pdbfile:
#            if residue.name in topmap.


if __name__ == '__main__':

    system = [('pope',162),('water',4860)]
    topology = {'pope': ['PEA','OLE','PAL']} # we're omitting water
    itp2top = {'pope': 'pope.itp'} 
    poten = PotentialEnergy('pope9x9chg_005.pdb',system,topology,itp2top)

# params = [sigma, epsilon]
# func_to_fit(params,configuration)
# fit_data = [(configuration,energy_value), (c1, ev1)]
# energy_value = lj+qq+lj14+qq14
# from Scientific.Functions import LeastSquares
# print LeastSqares.leastSqaresFit(func_to_fit,(.35,.3),fit_data)
