import Numeric as N

DEBYE = 4.8033324  # from electron-Ang to Debyes

def centerOfMass(coords,masses=None):
    if masses is None: masses = N.ones(len(coords))
    cms = N.add.reduce(coords*masses[:,N.NewAxis])/N.add.reduce(masses)
    return cms

def translateTo(coords,where):
    return coords-where

def dipoleMoment(coords,charges,masses=None):
    cms = centerOfMass(coords,masses)
    cms_coords = translateTo(coords,cms)
    dpl = N.add.reduce(cms_coords*charges[:,N.NewAxis])
    netdpl = N.sqrt(N.add.reduce(dpl*dpl))
    return netdpl*DEBYE
    
if __name__ == '__main__':
    # CH2
    coords = N.array([[13.340,7.440,46.080],
                      [14.360,7.800,45.920],
                      [12.760,8.330,46.320]])
    charges = N.array([-0.240,0.120,0.120])
    masses = N.array([12.011,1.008,1.008])
    print ' Dipole moment of CH2 is %.3g' % \
          dipoleMoment(coords,charges,masses)
    # CH3
    coords = N.array([[10.830, 12.460, 43.780],
                      [9.950, 12.230, 44.380],
                      [11.330, 13.240, 44.350],
                      [11.570, 11.660, 43.670]])
    charges = N.array([-0.480,0.160,0.160,0.160])
    charges = N.array([-0.180, 0.060, 0.06, 0.06])
    charges = N.array([-0.360, 0.120, 0.120, 0.120])
    masses = N.array([12.011,1.008,1.008,1.008])
    print ' Dipole moment of CH3 is %.3g' % \
          dipoleMoment(coords,charges,masses)
    # CH2-CH3
    coords = N.array([[10.320, 12.980, 42.440],
                      [10.830, 12.460, 43.780]])
    charges = N.array([-0.480,-0.240])
    masses = N.array([12.011,12.011])
    print ' Dipole moment of CH2-CH3 is %.3g' % \
          dipoleMoment(coords,charges,masses)
    # C-H
    coords = N.array([[10.830, 12.460, 43.780],
                      [9.950, 12.230, 44.380]])
    charges = N.array([-0.480,0.160])
    masses = N.array([12.011,1.008])
    print ' Dipole moment of C-H (terminal CH3 group) is %.3g' % \
          dipoleMoment(coords,charges,masses) 
    coords = N.array([[10.830, 12.460, 43.780],
                      [9.950, 12.230, 44.380]])
    charges = N.array([-0.240,0.120])
    masses = N.array([12.011,1.008])
    print ' Dipole moment of C-H (CH2 group) is %.3g' % \
          dipoleMoment(coords,charges,masses)
    coords = N.array([[10.830, 12.460, 43.780],
                      [9.950, 12.230, 44.380]])
    charges = N.array([-0.120,0.060])
    masses = N.array([12.011,1.008])
    print ' Dipole moment of C-H (CH2 group, original charges) is %.3g' % \
          dipoleMoment(coords,charges,masses)
    # C-C
    # to be done
    # CH3-CH2 total moment
    coords = N.array([[ 90.540,  18.190,  29.100],
                      [ 89.900,  17.330,  29.280],
                      [ 91.570,  17.840,  29.210],
                      [ 90.220,  19.250,  30.190],
                      [ 90.390,  18.770,  31.160],
                      [ 89.160,  19.490,  30.250],
                      [ 90.900,  20.090,  30.070]])
    #charges = N.array([-0.24, 0.12, 0.12, -0.48, 0.16, 0.16, 0.16])
    charges = N.array([-0.12, 0.06, 0.06, -0.18, 0.06, 0.06, 0.06])
    masses = N.array([12.011, 1.008, 1.008, 12.011, 1.008, 1.008, 1.008])
    print ' Dipole moment of CH2-CH3 is %.3g' % \
        dipoleMoment(coords,charges,masses)
    # CH2-CH2-CH3 
    coords = N.array([[ 90.370,  18.740,  27.630],
                      [ 91.020,  19.620,  27.600],
                      [ 89.350,  19.040,  27.390],
                      [ 90.540,  18.190,  29.100],
                      [ 89.900,  17.330,  29.280],
                      [ 91.570,  17.840,  29.210],
                      [ 90.220,  19.250,  30.190],
                      [ 90.390,  18.770,  31.160],
                      [ 89.160,  19.490,  30.250],
                      [ 90.900,  20.090,  30.070]])
    charges = N.array([ -0.24, 0.12, 0.12, 
                        -0.24, 0.12, 0.12, 
                        -0.48, 0.16, 0.16, 0.16])
    charges = N.array([ -0.12, 0.06, 0.06, 
                        -0.12, 0.06, 0.06, 
                        -0.18, 0.06, 0.06, 0.06])
    masses = N.array([12.011, 1.008, 1.008,
                      12.011, 1.008, 1.008, 
                      12.011, 1.008, 1.008, 1.008])
    print ' Dipole moment of CH2-CH2-CH3 is %.3g' % \
        dipoleMoment(coords,charges,masses)
    
