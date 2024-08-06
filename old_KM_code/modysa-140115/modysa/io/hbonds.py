import os.path, numpy as N
import modysa.misc as M
import cPickle, zipfile, ast, time

## http://mathesaurus.sourceforge.net/matlab-numpy.html octave/numpy

## FIXME: consider transforming HBTrajectory into more general
## *parent* class: EventTrajectory and use Inhertance to handle
## HBonds, WaterBridges, other interactions and events in time

class HBTrajectory:
    """store dictionaries {(donor,acceptor}: [(time_start,len),
    (time_start1, len1), ...]}; either ZipFile with entries named
    (time_block_start,time_block_end) and contents cPickle.dumps for
    dictionary data or map/dat file with zlib.compress'ed
    cPickle.dumps and stored (MAP) info about time_block_start, end
    and file position (tell) and byte count to be read; it seems
    reasonable to store HBonds for let's say 100 ps (block_size -
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
        self.step = step
        self.data = {}   # public data
        self.__data = {} # temporary data for fast browsing
        self.traj = None
        self.current = 0
        self.previous = None
        self.partno = -1

    def open(self, filename, mode='r'):
        assert mode in 'rwa'
        if mode is 'r':
            self.traj = zipfile.ZipFile(filename, mode)
            if hasattr(self.traj, 'comment'): 
                comment = ast.literal_eval(self.traj.comment)
                self.step = comment['Step']
            self.parts = [ast.literal_eval(x) for x in self.traj.namelist() \
                            if x[0] is '(' and x[-1] is ')']
            self.parts.sort()
        elif mode in 'wa':
            assert not os.path.lexists(filename)
            self.traj = zipfile.ZipFile(filename, mode, zipfile.ZIP_DEFLATED)
            comment = {'Created': time.ctime(), 'Step': self.step}
            self.traj.comment = str(comment)
            self.parts = []
            
    def close(self):
        assert self.traj
        # check if there are data not saved yet
        if self.traj.mode in 'wa' and self.data: self.sync()
        self.traj.close()

    def sync(self):
        assert self.traj and self.traj.mode in 'wa'
        if self.data:
            part_id = (self.previous, self.current)
            assert part_id not in self.parts
            # FIXME: check for overlaps! 
            # low, high = self.parts[idx]
            # self.previous >= low and self.current <= high
            self.traj.writestr(str(part_id),cPickle.dumps(self.data,2))
            #self.traj.fp.flush()
            #
            self.previous = self.current+self.step # first frame to be in the next part
        self.data = {}
        self.__data = {}

    def test(self):
        # self.parts continuous?
        # ravel, [1::2] - [:-1:2] = all_ones
        pass

    def add(self, data, frame_time):
        """ add new HBonds
        input: data is a list of (donor,acceptor) tuples
        HBonds are stored in self.data as 
        {(donor,acceptor): [single1, (start1, len1), single2..]"""
        position = int(frame_time/self.step)
        assert position > self.current # we *append* new HBonds
        if self.previous is None: self.previous = position
        self.current = position
        for pair in data:
            if self.data.has_key(pair):
                last_event = self.data[pair][-1]
                time_span = 0
                if isinstance(last_event, list):
                    first_seen, time_span = last_event
                else: first_seen = last_event
                if (position - first_seen - time_span) == 1:
                    time_span += 1
                    if time_span > 1: self.data[pair][-1][-1] = time_span
                    else: self.data[pair][-1] = [first_seen, time_span]
                else: self.data[pair].append(position)
            else: self.data[pair] = [position]

    # def periods(self):
    #     """convert raw data (self.data) from 
    #     [[start1,len1],start2,[start3,len3]] format into
    #     [[start1,len1],[start2,0],[start3,len3]]
    #     """
    #     if self.__data.has_key('periods'):
    #         if self.__data['periods']: return self.__data['periods']
    #     else: self.__data['periods'] = {}
    #     for pair, history in self.data.items():
    #         self.__data['periods'][pair] = self._unwrap(history)
    #     return self.__data['periods']

    def load(self, which=None):
        assert self.traj and 'r' in self.traj.mode
        if which is None: which = self.partno+1
        assert which<len(self.parts)
        self.data = cPickle.loads(self.traj.read(str(self.parts[which])))
        self.current = self.parts[which][0]
        self.partno = which
        # the contents of self.__data is needed to speed up browsing trajectory
        periods = {}
        for pair, history in self.data.items():
            periods[pair] = self._unwrap(history)
        # periods: [[40,2],[45,0],[80,6]]
        # events: [40,41,42,45,80,81,82,83,84,85,86]
        hb_id, history = zip(*periods.items())
        counts = N.array(\
            M.flatten([[i]*len(x) for i,x in enumerate(history)]), N.int)
        hits = N.array(M.flatten(history, full=False))
        hits[:,1] = hits[:,0] + hits[:,1]
        self.__data = {'counts': counts, 'hb_id': N.array(hb_id),
                       'hits': (hits[:,0], hits[:,1]), 'periods': periods}

    def history(self, donor, acceptor):
        # trajectory
        pair = (donor,acceptor)
        assert self.data and self.data.has_key(pair)
        history = self.data[pair]
        return self.step * N.array( M.flatten( \
                [isinstance(x,int) and [x] or \
                     range(x[0],x[0]+x[1]+1) for x in history]))
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

    def get(self, frame_time):
        """ return list of (donor, acceptor) tuples for a given frame_time 
        """
        assert 'r' in self.traj.mode # reading
        position = int(frame_time/self.step)
        idx = self.__locate(frame_time)
        if idx != self.partno: self.load(idx)
        #
        hit_low, hit_up = self.__data['hits']
        counts = self.__data['counts']
        hb_id = self.__data['hb_id']
        mask = N.logical_and(hit_low <= position, hit_up >= position)
        idx = counts[mask]
        hbonds = map(tuple,N.take(hb_id, idx, axis=0))
        return hbonds
        
        # in *load*
        # for pair, history in self.data.items():
        #     self.events[pair] = M.flatten([range(x[0],x[0]+x[1]+1) \
        #                         for x in self._unwrap(history)])
        ### events in format: [40,41,42,45,80,81,82,83,84,85,86]
        ### from self.data: [[40,2],45,[80,6]]
        # events = self.events
        # hbonds = [pair for pair in events.keys() \
        #               if position in events[pair]]

    def info(self):
        pass

    def limit(self,first=0,last=None):
        # FIXME - high priority!
        # load part containing first events
        # 1) _locate allowed idx range (first, last)
        # 2) trim marginal parts (lower, upper) - this migh be tricky
        ## modify get/filter (history?)
        # shrink self.data and self.__data in such a way
        # no earlier/later frames will be present in events
        pass

    # def count(self, histories):
    #     hb = [isinstance(x,int) and 1 or x[-1]+1 for x in histories]
    #     ## for x in M.flatten(histories)]
    #     return sum(hb)
    def count(self, events):
        if len(events) is 0: return 0
        return N.add.reduce(N.array(events)[:,1]+1)
            
    def filter(self, subset, acceptor=False):
        """ return list of events (two-element lists) 
        involving atoms from a subset """
        assert self.data
        assert acceptor in [True, False]
        # FIXME: some inconsistency here (cf. get)
        # returned *flatten* list of events (i.e. [[start,len],start])
        # or list of periods (i.e. [[start,end],[start,end]]
        # maybe we could use self.__data here to avoid some array operations?
        unique = N.unique1d(subset) # donors by default
        hb_id, hb_history = zip(*self.data.items())
        atoms = zip(*hb_id)[acceptor]
        # http://stackoverflow.com/questions/1273041/how-can-i-implement-matlabs-ismember-command-in-python
        auniq, ainv = N.unique1d(atoms, return_inverse=True)
        idx = N.where(N.take(N.setmember1d(auniq, unique), ainv))[0]
        nd_history = N.array(hb_history,N.object)
        return self._unwrap(M.flatten(N.take(nd_history, idx).tolist()))
        # see also:
        # http://stackoverflow.com/questions/1613249/numpy-comparing-elements-in-two-arrays

