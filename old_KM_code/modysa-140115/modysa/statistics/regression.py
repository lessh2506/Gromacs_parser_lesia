import numpy as N

class LinearRegression:

    def __init__(self, x, y):
        """calculate linear regression: y = b0 + b1*x  """

        n = len(x)
        sumx = N.add.reduce(x)
        sumy = N.add.reduce(y)
        sumxx = N.add.reduce(x*x)
        sumxy = N.add.reduce(x*y)
        sumyy = N.add.reduce(y*y)
        detC = n*sumxx - sumx*sumx
        self.b0 = (sumy*sumxx - sumxy*sumx)/detC    # intercept
        self.b1 = (n*sumxy - sumx*sumy)/detC        # slope
        varx = (sumxx - sumx*sumx/n)/(n-1)            # variance of x
        vary = (sumyy - sumy*sumy/n)/(n-1)            # variance of y
        avex = sumx/n
        avey = sumy/n
        self.r = self.b1*N.sqrt(varx)/N.sqrt(vary)           # correlation coef
        varr = (n-1)*(vary - self.b1*self.b1*varx)/(n-2)      # residual variance
        self.sb0 = N.sqrt(varr*sumxx/detC)            # intercept stddev
        self.sb1 = N.sqrt(varr*n/detC)                # slope stddev
        tb0 = self.b0/self.sb0                                # t-values
        tb1 = self.b1/self.sb1

    def __str__(self):
        # print nicely results
        out = [' results of fitting data to y = b0 + b1 * x ',
               ' b0 (sb0): %12.5g (%12.5g)' % (self.b0, self.sb0),
               ' b1 (sb1): %12.5g (%12.5g)' % (self.b1, self.sb1),
               ' R       : %12.5g' % self.r]
        return '\n'.join(out)

    def __call__(self):
        return self.b0, self.b1, self.sb0, self.sb1, self.r
