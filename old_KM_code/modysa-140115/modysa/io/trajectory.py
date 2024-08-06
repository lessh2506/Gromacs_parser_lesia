import glob, re, numpy as N
import modysa.misc as M
from gromacs import xtcio as X



class Trajectory:

    __errors = {'known_type': 'Trajectory type not yet supported',
                'unknown_type': 'Trajectory type not supported '
                                'or bad filename pattern'}

    def __init__(self, filename_pattern):

        self.filenames = glob.glob(filename_pattern)
        if not self.filenames: raise IOError, 'No trajectory files found!'
        self.filenames.sort()

        self.current = 0
        self.parts = []

        self.dt = None         # time between xyz dumps
        self.time_step = None  # MD delta_t, usually 2 fs
        self.last = None
        self.step = None

        t_type = self.type()
        if t_type == 'XTC':
            self.__handleXTC()
        elif t_type == 'TRR':
            raise ValueError(self.__errors['known_type'])
        elif t_type == 'XYZ': # amber
            raise ValueError(self.__errors['known_type'])
        else:
            raise ValueError(self.__errors['unknown_type'])

        # finishing touches
        self.reset()

    def __handleXTC(self): 
        self.traj = X.XTCReader()
        # find trajectory sampling
        self.traj.open(self.filenames[0])
        self.traj.readFrame()
        t0 = self.getTime()
        s0 = self.getStep()
        self.traj.readFrame()
        t1 = self.getTime()
        s1 = self.getStep()
        self.step = s1 - s0
        self.dt = t1 - t0
        self.time_step = self.dt/float(self.step)

        for fnm in self.filenames: # [1:]?
            self.traj.open(fnm)
            t0 = self.traj.getTime()
            tn = self.traj.getLastFrameTime()
            self.traj.close()
            self.parts.append((t0, tn))
        # check for overlaps in consecutive parts
        # if there are overlaps, the contents of i+1 overrides the ith part
        for i in range(len(self.parts)-1):
            self.parts[i] = (self.parts[i][0], self.parts[i+1][0]-self.dt)
        # the last part remain untouched: self.part[-1] = self.part[-1]

    def reset(self):
        self.first = self.parts[0][0]
        self.last = self.parts[-1][-1]
        self.range = range(len(self.parts))
        self.current = self.range[0]
        self.now = self.first
        self.open()

    def type(self):
        assert self.filenames
        ext_pattern = re.compile('^.+\.(.+)$', re.M)
        filenames = '\n'.join(self.filenames)
        hits = re.findall(ext_pattern, filenames)
        if len(self.filenames) == len(hits): 
            # OK, all filenames have the same extension
            return hits[0].upper()

    def limit(self, first=None, last=None):
        assert self.dt is not None

        if first is None and self.first is not None: first = self.first
        if last is None and self.last is not None: last = self.last
        #last += self.dt
        assert first < last
        self.range = [i for i in range(len(self.parts)) \
                          if self.parts[i][1] >= first and \
                             self.parts[i][0] < last]
        assert self.range
        self.first = first
        self.last = last
        self.current = self.range[0]
        self.open()
        self.traj.seekByTime(self.first)
        self.now = self.getTime()

    def open(self):
        assert self.traj
        # FIXME: open if not opened previously
        self.traj.open(self.filenames[self.current])
        self.now = self.getTime()

    def close(self):
        assert self.traj
        if self.traj._XTCReader__fileDescr != -1:
            self.traj.close()

    def getTime(self):
        assert self.traj
        t = self.traj.getTime()
        if t: return M.floatSignificantNumbers(t)
        else: return 0 # None ...

    def getStep(self):
        assert self.traj
        return self.traj.getStep()

    def getCell(self):
        assert self.traj
        return self.traj.getCell()

    def getConfiguration(self):
        assert self.traj
        return self.traj.getCoords()

    def readFrame(self):
        assert self.traj


        def checkGap(prev):
            self.now = self.getTime()
            if (self.now-prev) > self.dt:
                raise ValueError('A gap in trajectory detected (%d, %f, %f)' \
                                     % (self.current, self.now, prev))

        if self.dt is None and len(self.parts) <= self.current:         
            # on opening trajectory no self.dt yet
            self.traj.readFrame()
            return True
        else:
            t = self.getTime()
            if t < self.parts[self.current][-1]: 
                self.traj.readFrame()
                checkGap(t)
            else:
                self.current += 1
                self.close()
                if self.current in self.range:
                    # opening implies reading the first frame
                    self.open()
                    checkGap(t)
                else: return False # end of trajectory reached
            if self.last is not None and self.now <= self.last: 
                return True # ready to read next
            else: return False # end of set limit reached

