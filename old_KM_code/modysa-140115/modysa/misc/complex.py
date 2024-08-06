
import numpy as np
import cmath

class PeriodicDataSeries:
    # values in radians, i.e. in (-pi/pi) range
    # learnt from Q Rev Biophys (1998) 31(2):145-237
    # Peter Guntert, Structure calculation of biological macromolecules 
    # from NMR data, p. 7.3, pp. 205
    # ave = \arg \sig_k \exp (j \phi_k)
    # stdev = -2. \log | \frac{1}{k} \sig_k \exp (j \phi_k)
    def __init__(self, data):
        self.data = np.array(data)
        
    def mean(self):
        return cmath.phase(np.sum(np.exp(1j*self.data)))

    def standardDeviation(self):
        ave = np.sum(np.exp(1j*self.data))/float(len(self.data))
        magnitude = np.sqrt(ave.real**2 + ave.imag**2)
        return -2.*np.log(magnitude)

rad2deg = lambda value: value*180./np.pi
deg2rad = lambda value: value*np.pi/180.

if __name__ == '__main__':
    sample = [-175., 174, 178, -170, 179]
    a = PeriodicDataSeries(map(deg2rad,sample))
    print 'ave: %12.3f stdev: %12.3f' % tuple(map(rad2deg,
                                 (a.mean(), a.standardDeviation())))
