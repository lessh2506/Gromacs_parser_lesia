
import numpy as N
# Numerical 23.8-1 bug: len(n[1:]) != len(n[:-1])
import cmath
import modysa.misc as M
import os.path as OP, inspect
import scipy.optimize as SOPT


class Histogram:

    def __init__(self, periodic=False):
        assert periodic in [False,True]
        self.__periodic = periodic
        self.bins = None

    def load(self, filename):
        # FIXME: read bins/data
        bins, counts = N.loadtxt(filename).transpose()
        self.bins = bins
        self.n = len(self.bins)
        self.data = counts

    def save(self, filename, fmt=None):
        # FIXME: os.path.lexists
        assert not OP.lexists(filename), 'File %s exists. Rename.' % filename
        if fmt is None: fmt = '%15.5e %15.5e'
        N.savetxt(filename, zip(self.bins, self.data), fmt=fmt)

    def set(self, low=None, up=None, nbins=None, bins=None):
        assert None not in [low,up,nbins] or bins is not None
        if bins is not None: self.bins = N.array(bins)
        else: self.bins = N.linspace(low, up, nbins)
        self.n = len(self.bins)
        if self.__periodic and self.bins[0] < -N.pi and self.bins[-1] > N.pi:
	    self.__binscale = N.pi/180.
        else: self.__binscale = 1.0
        self.reset()

    def reset(self):
        self.data = N.zeros(len(self.bins), N.ulonglong)

    def update(self, other):
        assert isinstance(other, Histogram)
        assert other.n == self.n
        self.data += other.data

    def addData(self, data, weights=None):
        # FIXME: if weights are not None: each data point has its weight
        # see Scientific.Statistics
        n = N.searchsorted(N.sort(data),self.bins)
        n = N.concatenate([n,[len(data)]])
        self.data += n[1:] - n[:-1]

    def __call__(self, norm=1):
        count = N.add.reduce(self.data)
        if norm != 0: factor = float(norm)/count
        else: factor = 1.
        # get some statistics
        # *** mean ***
        if not self.__periodic: 
            ave = sum(1.0*self.bins*self.data)/sum(self.data)
        else: 
            nz = N.nonzero(self.data)[0]
            # avoid operations on bins with 0 weight
            bins = self.bins[nz] * self.__binscale
            res = N.sum(self.data[nz]*N.exp(1j*bins)) / \
                N.sum(self.data[nz])
            ave = N.angle(res)/self.__binscale
        self.mean = ave
        # *** standardDeviation ***
        if not self.__periodic:
            mean2 = N.sum(1.0 * self.bins**2 * self.data)/sum(self.data)
            std = N.sqrt(mean2 - self.mean**2)
        else:
            std = N.sqrt(1-abs(res))/self.__binscale
        self.standardDeviation = std
        return factor*self.data

    def samples(self, total=None):
        # FIXME: if histogram is normalized use total to regenerate counts
        nz = N.nonzero(self.data)[0]
        s = M.flatten([[i]*c for i, c in zip(self.bins[nz],self.data[nz])])
        return s


class Fit:

    # http://www.scipy.org/Cookbook/FittingData
    fitfunc = None
    x = None

    def __init__(self):
        self.fitted = None # fitted data
        self.params = None
        self.errfunc = lambda params, x, data: self.fitfunc(params, x) - data

    def __call__(self, data, *p0):
        assert None not in [self.fitfunc, self.x]
        if not p0: 
            # FIXME: consider .defaults .keywords
            nargs = len(inspect.getargspec(self.fitfunc).args)
            # find the number of arguments for the fitted function, set'em to one
            p0 = N.ones(nargs)
        p1, success = SOPT.leastsq( self.errfunc, p0[:], args=(self.x, data))
        self.params = p1
        self.fitted = self.fitfunc(p1, self.x)

    def stats(self):
        # report quality of the fit
        # mean, stdev, mode
        pass


class FitPeriodicNormalDistribution(Fit):

    # assume von Mises distribution (wrapped normal dist)
    fitfunc = None


class FitNormalDistribution(Fit):
    fitfunc = None


class myFit(Fit):

    @staticmethod
    def fitfunc(params, x):
        return N.sin(x) * params[0] * N.exp(-.5*(x/params[1])**2)

    x = N.linspace(0, N.pi, 1800) # fit only in (0,60) deg range

def mean(data):
    return N.mean(data)

def standardDeviation(data):
    return N.std(data, ddof=1)

# stolen from http://www.scipy.org/Cookbook/SignalSmooth
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; 
            should be an odd integer
        window: the type of window from 'flat', 'hanning', 
            'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself 
          if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', \
                          'hanning', 'hamming', 'bartlett', 'blackman'"


    s=N.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=N.ones(window_len,'d')
    else:
        w=eval('N.'+window+'(window_len)')

    y=N.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

class DataSeries:

    def __init__(self,data):
        self.data = N.array(data)
        self.reset()

    def reset(self):
        for item in ['mean','rho','f','n',
                     'standardDeviation','standardError']:
            setattr(self, item, None)

    def _mean(self, data):
        return N.mean(data,axis=-1)

    def _standardDeviation(self, data):
        return N.std(data,axis=-1)

    def __call__(self, block_size=None):
        if block_size: 
            block_size = int(block_size)
            assert block_size < len(self.data)
            nblocks = len(self.data)/block_size
            offset = len(self.data) - nblocks*block_size
            data = self._mean(self.data[offset:].reshape(-1,block_size))
            data = N.array(data)
        else: data = self.data

        self.mean = self._mean(data)
        self.n = len(data)
        self.standardDeviation = self._standardDeviation(data)

        if self.n > 2:
            series = data - self.mean
            self.rho = sum(series[1:]*series[:-1])/sum(series[1:-1]**2)
            if self.rho > 1: self.rho = 0.99   # see bence, 1995
            elif self.rho < -1: self.rho = -.99
            f = (1.+self.rho)/(1.-self.rho)
            if N.isnan(f): self.f = 1.
            else: self.f = N.sqrt(f)
        else:
            self.f = 1.
            self.rho = None

        self.standardError    = self.f*self.standardDeviation/N.sqrt(self.n)

class PeriodicDataSeries(DataSeries):
    # values in radians, i.e. in (-pi/pi) range
    # learnt from Q Rev Biophys (1998) 31(2):145-237
    # Peter Guntert, Structure calculation of biological macromolecules 
    # from NMR data, p. 7.3, pp. 205
    # ave = \arg \sig_k \exp (j \phi_k)
    # stdev = -2. \log | \frac{1}{k} \sig_k \exp (j \phi_k)

    def _mean(self, data):
        data = N.sum(N.exp(1j*N.array(data)),axis=-1)
        if M.isiterable(data): mean = map(cmath.phase, data)
        else: mean = cmath.phase(data)
        return mean

    def _standardDeviation(self, data):
        # FIXME: sanity check needed, array type/shape
        imean = N.sum(N.exp(1j*N.array(data)))/float(len(data))
        magnitude = N.sqrt(imean.real**2 + imean.imag**2)
        return -2.*N.log(magnitude)

