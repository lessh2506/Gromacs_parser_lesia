import numpy as N


EMPTY_TAG = '-'

class StatePath:

    def __init__(self, path=None):
        self.set(path)

    def __update(self, other):
        path = [x is EMPTY_TAG and y or x for x,y in zip(self.path,other)]
        return ''.join(path)
        
    def __add__(self, other):
        return StatePath(self.__update(other))

    def __str__(self):
        return self.path

    def set(self, path):
        if path is None: path = ''
        else: assert type(path) is type('')
        self.path = path

    def update(self, other):
        self.path = self.__update(other)

class DihedralAnalysis:

    named_conformations = {\
        'trans': lambda a: N.logical_or(a < -150, a > 150),
        'gauche_minus': lambda a: N.logical_and(a > -90, a < -30),
        'gauche_plus': lambda a: N.logical_and(a > 30, a < 90),
        'gauche': lambda a: N.logical_and(N.abs(a) > 30, N.abs(a) < 90)
        }

    named_tags = {'trans': 't', 'gauche_minus': 'm',
                  'gauche_plus': 'p', 'gauche': 'g'}

    def getConformations(self, angle, 
                         limits=None, named=None, 
                         mode=None, tag='x'):
        """ return list of values in a given conformation,
        you can either give the name of a particular conformation
        (see named_conformations) or limit dihedral values
        defining a conformation; in the latter case you can set limits with
        either two-element-lists (angle < low or angle > high) or
        two-element-tuple (low < angle < high)
        """
        if type(angle) in [list, tuple]: angle = N.array(angle)
        assert (limits is not None and named is None) or \
            (limits is None and named is not None), \
            'Pick either _limits_ or _named_'
        modes = ['which', 'values', 'path']
        if mode is None: mode = modes[0]
        assert mode in modes
        if named:
            assert self.named_conformations.has_key(named)
            func = self.named_conformations[named]
        else:
            assert len(limits) == 2
            low, high = limits
            if type(limits) == list: # FIXME: mutable.. 
                func = lambda a: N.logical_or(a < low, a >= high)
            elif type(limits) == tuple:
                func = lambda a: N.logical_and(a >= low, a < high)
            else: func = lambda a: N.ones(len(a), N.int)
        idx = N.flatnonzero(func(angle))
        result = None
        if mode is modes[0]:   result = idx.tolist()
        elif mode is modes[1]: result = N.take(angle, idx).tolist()
        elif mode is modes[2]: 
            result = N.chararray(len(angle),buffer=EMPTY_TAG*len(angle))
            if named: tag = self.named_tags[named]
            N.put(result, idx, tag)
            result = ''.join(result)
        return result

    def mean(self, angle):
        # imaginary
        pass

    def standardDeviation(self, angle):
        pass
