from Scientific.IO import NetCDF
import os.path as OP, numpy as N, marshal, binascii, cPickle, zipfile
import modysa.misc as M

class AngularTrajectory:

    """ the idea taken from nmoldyn.core.AngularTrajectory and 
        nmoldyn.calc.qTrajectory"""

    def __init__(self):
        # self.step = step
        self.__reset()

    def __reset(self):
        for item in ['traj', 'mode', 'current', 'description', 'valid']:
            setattr(self, item, None)

    def open(self, filename, mode='r'):
        assert mode in 'rw'
        if mode is 'w': assert not OP.lexists(filename) or overwrite
        self.traj = NetCDF.NetCDFFile(filename, mode)
        self.mode = mode
        self.current = 0
        if mode is 'r':
            self.valid = self._isvalid()
            self.description = \
                marshal.loads(binascii.a2b_base64(self.traj.description))
            # description should contain definition of atomic groups
        elif mode is 'w':
            # nmoldyn.core.AngularTrajectory
            self.traj.file.createDimension('quaternion_length',4)
            pass
        self.valid = False

    def close(self):
        self.traj.close()

    def write(self, data, frame_time=None):
        # MMTK.Trajectory.RigidBodyTrajectory
        pass

    def read(self, item, first=0, last=None, skip=1, varname='com'):
        # varname: com (center-of-mass), quaternions, mass, rmsd (fit quality)
        # note that one can access these variables directly:
        # self.traj.file.variables[varname][item,first:last:skip]
        pass



# see density.DensityTrajectory

class DataTrajectory:

    def __init__(self, dump_step_time=1.0):
        self.__reset()
        self.dump_step_time = dump_step_time

    def __reset(self):
        for item in ['traj', 'mode', 'current', 'names', 'valid']:
            setattr(self, item, None)

    def read(self, *kinds, **limits):
        start, end, step, flat = 0, None, None, False
        for l,v in limits.items():
            if l == 'start': start = v
            elif l == 'end': end = v
            elif l == 'step': step = v
            elif l == 'flat': flat = v
        assert self.traj is not None
        result = None
        for k in kinds:
            assert self.names.has_key(k)
            gj = self.traj.variables['data'][start:end:step,self.names[k]]
            if result is None: result = gj
            else: result += gj
        if flat: # FIXME: mean or sum?
            result = N.add.reduce(result)/len(result)
        return result

    def open(self, filename, mode='r', overwrite=False):
        assert mode in 'rw'
        if mode is 'w': assert not OP.lexists(filename) or overwrite
        self.traj = NetCDF.NetCDFFile(filename, mode)
        self.mode = mode
        self.current = 0
        if mode is 'r':
            self.valid = self._isvalid()
            self.names = marshal.loads(binascii.a2b_base64(self.traj.names))
        self.valid = False

    def _isvalid(self):
        if self.traj is None or not hasattr(self.traj, 'dimensions') or \
            not hasattr(self.traj, 'variables') or \
            not hasattr(self.traj, 'names'): return False
        for dim in ['size', 'kind', 'step']:
            if not self.traj.dimensions.has_key(dim): return False
        for var in ['data', 'time']:
            if not self.traj.variables.has_key(var): return False
        if not hasattr(self.traj, 'names'): return False
        return True

    def close(self):
        self.traj.close()

    def write(self, data, frame_time=None, scale=True):
        if not self.valid:
            labels = sorted(data.keys())
            assert self.traj, 'No valid DensityTrajectory opened!'
            self.traj.createDimension('step', None)
            # FIXME: what if various items have different size?
            # no kind/size: instead ranges of indices, e.g.
            # for 12-C, 14-C, 10-C, 14-C chains:
            # (0,12), (12,26), (26,36), (36,50)
            # thus in names we should store ranges not ordinal numbers and
            # there should be only two dimensions (step and size, with size
            # being sum of all items (50 in an example above, range: <0-49>)
            self.traj.createDimension('size', len(data.items()[0][1]))
            self.traj.createDimension('kind', len(data))
            self.names = dict([(labels[i], i) for i in range(len(labels))])
            self.traj.names = binascii.b2a_base64(marshal.dumps(self.names))
            # self.traj.history = "Created ..."
            # self.traj.description = "Density profile collected.."
            self.traj.createVariable('data', 'f', ('step', 'kind', 'size'))
            self.traj.createVariable('time', 'd', ('step',))
            self.valid = self._isvalid()

        if frame_time is None: frame_time = self.current * self.dump_step_time
        self.traj.variables['time'][self.current] = frame_time

        assert len(data) == self.traj.dimensions['kind']
        for k in data.keys():
            if type(data[k]) != type(N.array([1])): d = N.array(data[k])
            else:                                   d = data[k]
            if scale and type(scale) == int:        d *= scale
            self.traj.variables['data'][self.current, self.names[k]] = \
                d.astype('f')
        self.traj.sync()
        self.current += 1

## http://mathesaurus.sourceforge.net/matlab-numpy.html octave/numpy

class EventTrajectory:
    """container for various *events* (HBonds, AtomContacts, WaterBridges); 
    events are identified as tuples of atom indices (two-, three-element, etc.)

    store dictionaries {(atomA,atomB}: [(time_start,len),
    (time_start1, len1), ...]}; either ZipFile with entries named
    (time_block_start,time_block_end) and contents cPickle.dumps for
    dictionary data or map/dat file with zlib.compress'ed
    cPickle.dumps and stored (MAP) info about time_block_start, end
    and file position (tell) and byte count to be read; it seems
    reasonable to store *events* for let's say 100 ps (block_size -
    parameter defined by a user) since after loading such a chunk (not
    too small since performance will be hurted by disk IO and not too
    large since it has to fit into memory) we might process data using
    efficient numpy calls;

    e.g. let's say we have a list_of_donors, list_of_acceptors and
    list_of_HBhistory, let's say we also have allowed_donors and
    allowed_acceptors, then we can make

    d = np.where(np.setmember1d(list_of_donors,allowed_donors))[0]
    a = np.where(np.setmember1d(list_of_donors,allowed_acceptors))[0]
    hb = np.where(np.setmember1d(d,a))[0]
    allowed_HBonds = np.take(list_of_HBhistory,np.take(d,hb)) 

    d and a are indices, hb as well; allowed_HBonds are list of HB histories

    it's reasonable to store HBond data as dictionary with (donor,acceptor)
    keys since we can easily extract particular HBond and if we want
    to process many HBonds at the same time (the numpy recipe rescribed 
    above) we might easily convert dictionary into three lists:

    hb_id, list_of_HBhistory = data.items()
    list_of_donors, list_of_acceptors = zip(*hb_id)"""

    def __init__(self, step=1.0):
        # FIXME: step to open('wa')
        # FIXME: step read from info while open('r')
        self.step = float(step)
        self.current = 0
        self.previous = None

        self.data = {}   # public data
        self.__data = {} # temporary data for fast browsing
        self.traj = None
        self.partno = -1

    def open(self, filename, mode='r'):
        assert mode in 'rwa'
        if mode is 'r':
            self.traj = zipfile.ZipFile(filename, mode)
            self.parts = [eval(x) for x in self.traj.namelist() \
                            if x[0] is '(' and x[-1] is ')']
            self.parts.sort()
            self.current = self.parts[0][0]
        elif mode in 'wa':
            assert not OP.lexists(filename)
            self.traj = zipfile.ZipFile(filename, mode, zipfile.ZIP_DEFLATED)
            self.parts = []
            
    def close(self):
        assert self.traj
        # check if there are data not saved yet
        if self.traj.mode in 'wa' and self.data: self.sync()
        self.traj.close()

    def sync(self):
        assert self.traj and self.traj.mode in 'wa'
        if self.data:
            part_id = (int(self.previous), int(self.current))
            assert part_id not in self.parts # has_key?
            # FIXME: check for overlaps! 
            # low, high = self.parts[idx]
            # self.previous >= low and self.current <= high
            self.traj.writestr(str(part_id),cPickle.dumps(self.data,2))
            #self.traj.fp.flush()
            #
            # first frame to be in the next part
            self.previous = self.current + 1
        self.data = {}
        self.__data = {}

    def test(self):
        # self.parts continuous?
        # ravel, [1::2] - [:-1:2] = all_ones
        pass

    def add(self, data, frame_time):
        """ add new events
        input: data is a list of (atomA,..,atomN) tuples, 
        e.g. HBonds are stored in self.data as 
        {(donor,acceptor): [single1, (start1, len1), single2..]"""
        position = int(frame_time/self.step)
        assert position > self.current # we *append* new events
        if self.previous is None: self.previous = position
        self.current = position
        for atuple in data:
            if self.data.has_key(atuple):
                last_event = self.data[atuple][-1]
                time_span = 0
                if isinstance(last_event, list):
                    first_seen, time_span = last_event
                else: first_seen = last_event
                if (position - first_seen - time_span) == 1:
                    time_span += 1
                    if time_span > 1: self.data[atuple][-1][-1] = time_span
                    else: self.data[atuple][-1] = [first_seen, time_span]
                else: self.data[atuple].append(position)
            else: self.data[atuple] = [position]

    def __len__(self):
        assert self.traj and 'r' in self.traj.mode
        return self.parts[-1][1]-self.parts[0][0] + 1

    def load(self, which=None):
        assert self.traj and 'r' in self.traj.mode
        if which is None: which = self.partno+1
        assert which<len(self.parts)
        self.data = cPickle.loads(self.traj.read(str(self.parts[which])))
        self.current = self.parts[which][0]
        self.partno = which
        # the contents of self.__data is needed to speed up browsing trajectory
        periods = {}
        for atuple, history in self.data.items():
            periods[atuple] = self._unwrap(history)
        # periods: [[40,2],[45,0],[80,6]]
        # events: [40,41,42,45,80,81,82,83,84,85,86]
        event_ids, history = zip(*periods.items())
        counts = N.array(\
            M.flatten([[i]*len(x) for i,x in enumerate(history)]), N.int)
        hits = N.array(M.flatten(history, full=False))
        hits[:,1] = hits[:,0] + hits[:,1]
        self.__data = {'counts': counts, 'event_ids': N.array(event_ids),
                       'hits': (hits[:,0], hits[:,1]), 'periods': periods}

    def history(self, *atuple):
        # PREVIOUS: history(self. donor, acceptor)
        # trajectory
        assert self.data and self.data.has_key(atuple)
        history = self.data[atuple]
        return M.flatten([isinstance(x,int) and [x] or \
                              range(x[0],x[0]+x[1]+1) for x in history])
        #return self._history2events(self.data[pair])

    def __locate(self, frame_time):
        position = int(frame_time/self.step)
        # FIXME: store self.parts as zip(*list_of_tuples)
        nd_parts = N.array(self.parts)
        idx = N.searchsorted(nd_parts[:,1], position)
        assert idx < len(self.parts), \
            'frame_time: %g, not found (last time range: %s)' % \
            (frame_time, str((self.step*self.parts[-1][0],
                              self.step*self.parts[-1][1])))
        # FIXME: the condition below is not necessery?
        # assert position >= self.parts[idx][0] \
        #     and position <= self.parts[idx][1], \
        #     'frame_time: %g, not found (tried %s time range)' % \
        #     (frame_time, str((self.step*self.parts[idx][0],
        #                       self.step*self.parts[idx][1])))
        return idx

    def _unwrap(self,history):
        # gj = []
        # for history in self.data.values():
        #     gj.extend([isinstance(x,int) and 1 or 0 for x in history])
        ## sum(gj) -> 20970, len(gj) -> 43296
        ## should we use more memory consuming but faster format 
        ## [[start1,len1], [start2, len2]...]

        return [isinstance(x,int) and [x,0] or x for x in history]

        #return M.flatten([isinstance(x,int) and [x] or \
        #                            range(x[0],x[0]+x[1]+1) for x in history])

    def get(self, frame_time=None):
        """ return list of tuples of atom indices for a given frame_time 
        """
        assert 'r' in self.traj.mode # reading
        if frame_time is None: frame_time = self.current * self.step
        if frame_time > self.parts[-1][-1]: return None
        position = int(frame_time/self.step)
        idx = self.__locate(frame_time)
        if idx != self.partno: self.load(idx)
        #
        hit_low, hit_up = self.__data['hits']
        counts = self.__data['counts']
        event_ids = self.__data['event_ids']
        mask = N.logical_and(hit_low <= position, hit_up >= position)
        idx = counts[mask]
        events = map(tuple,N.take(event_ids, idx, axis=0))
        self.current = position
        return events
    
    def next(self):
        self.current += 1
        return self.get()
    
    def info(self):
        pass

    def time(self):
        # frame_time
        return float(self.current) * self.step

    def limit(self,first=0,last=None):
        # FIXME - high priority!
        # load part containing first events
        # 1) _locate allowed idx range (first, last)
        # 2) trim marginal parts (lower, upper) - this migh be tricky
        ## modify get/filter (history?)
        # shrink self.data and self.__data in such a way
        # no earlier/later frames will be present in events
        pass

    def count(self, events):
        if len(events) is 0: return 0
        return N.add.reduce(N.array(events)[:,1]+1)
            
    def filter(self, subset, which=0):
        """ return list of event_ids (tuples of atom indices)
        involving atoms from a given subset (*which* defines the position
        of the atom index in a tuple)"""
        if not self.data and 'r' in self.traj.mode: self.load()
        else: assert self.data

        # FIXME: some inconsistency here (cf. get)
        # returned *flatten* list of events (i.e. [[start,len],start])
        # or list of periods (i.e. [[start,end],[start,end]]
        # maybe we could use self.__data here to avoid some array operations?
        unique = N.unique1d(subset) 
        event_ids, history = zip(*self.data.items())
        atoms = zip(*event_ids)[which]
        # http://stackoverflow.com/questions/1273041/how-can-i-implement-matlabs-ismember-command-in-python
        auniq, ainv = N.unique1d(atoms, return_inverse=True)
        idx = N.where(N.take(N.setmember1d(auniq, unique), ainv))[0]
        nd_history = N.array(history,N.object)
        return self._unwrap(M.flatten(N.take(nd_history, idx).tolist(),full=False))

        # see also:
        # http://stackoverflow.com/questions/1613249/numpy-comparing-elements-in-two-arrays

