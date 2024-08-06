
import numpy as N, string
import itertools, operator, re


class Units:
    deg = N.pi/180. # rad
    cm  = 1e-2      # m
    nm  = 1e-9      # m
    cal = 4.184     # kJ
    
class Constants:
    N_avo = 6.022141793e23  # Avogadro constant

# def flatten(seq):
#     return tuple(itertools.chain.from_iterable(seq))

def flatten(seq, full=True, maxdepth=10):
    done = True
    idepth = 0
    while done:
        seq = list(itertools.chain.from_iterable(\
                [ isiterable(x) and x or [x] for x in seq]))
        done = not isflat(seq) and full
        assert idepth < maxdepth, 'Error: maximum number of nesting reached'
        idepth += 1
    return seq
        
def setDefault(params, name, default_value):
    if params.has_key(name): return params[name]
    else: return default_value

def comb2el(set):
    # python2.6: itertools.combinations(set,2)
    a = set[:-1]
    b = set[1:]
    for i in range(2,len(set)):
        a.extend(set[:-i])
        b.extend(set[i:])
    return zip(a,b)

class Results(dict):

    def update(self, data):
        for k,v in data: 
            if self.has_key(k): 
                if type(v) == type([]): self[k].extend(v) # iterable
                else: self[k].append(v)
            else: 
                if type(v) == type([]): self[k] = v
                else: self[k] = [v]

    def setUp(self, keys):
        for k in keys: self[k] = []

    def reset(self):
        for k in self.keys(): self[k] = []

def occurence(data):
    ## taken from: http://stackoverflow.com/questions/4651683/numpy-grouping-using-itertools-groupby-performance
    # values = N.array(data)
    # values.sort()
    # diff = N.concatenate(([1], N.diff(values)))
    # idx = N.concatenate((N.where(diff)[0], [len(values)]))
    # index = N.empty(len(idx)-1)
    # return zip(values[idx[:-1]], N.diff(idx))
    ## http://www.builderau.com.au/program/python/soa/Python-groupby-the-iterator-swiss-army-knife/0,2000064084,339280431,00.htm
    ## sum(1 for i in v) is supposed to be 2x faster than len(list(v))
    return dict([(k,sum(1 for i in v)) for k,v in \
                     itertools.groupby(sorted(data))])

def slicer(data):
    # identify consequtive runs of integers
    # e.g. from [1,  4,5,6,  10,  15,16] get [1], [4,5,6], [10], [15,16]
    # FIXME: add parameter defining step between neighbor integers
    #        (default: 1)
    return [map(operator.itemgetter(1), g) \
                  for k, g in itertools.groupby( \
            enumerate(sorted(data)), lambda (index, item): index-item)]
    
# pick_items function takes limiting values from list of (start, end)
# and returns elements of data array found inside of such intervals
# beware of the fact that data *must* be a numpy array not a list
# EXAMPLE: find distances from intervals [(2.,4.5), (6.,8.)]
#          pick_items( distances, intervals )

def pick_items(data, intervals):
    assert isiterable(data) and isiterable(intervals)
    orig_type = type(data)
    if orig_type != N.ndarray: data = N.array(data)
    result = N.nonzero(\
        reduce(N.logical_or, [N.logical_and(data>=a_min,data<a_max) \
                                  for a_min, a_max in intervals]))[0]
    if orig_type != N.ndarray: return orig_type(result)
    else: return result # N.ndarray(shape=result.shape, buffer=result, dtype=result.dtype)

#pick_items = lambda data, intervals: N.nonzero(\
#   reduce(N.logical_or, [N.logical_and(data>=a_min,data<a_max) \
#                             for a_min, a_max in intervals]))

# generate indices from a list of intervals <start,end)
mrange = lambda intervals: flatten(map(lambda x: range(*x), intervals))


def isiterable(seq, excludedtypes=(str,)):
    # default value of excludedtypes allow to call 
    # flatten function on nested lists of strings
    try: 
        result = iter(seq) and not type(seq) in excludedtypes and True
    except TypeError: result = False
    return result

def isflat(seq):
    x = sum(map(isiterable,seq))
    if x == 0: return True
    else: return False

def iscallable(instance):
    return hasattr(instance, '__call__')
    

def parseStr(x):
    # convert a string into a int or float or if none of those return a string
    return x.isalpha() and x or \
           x.isdigit() and int(x) or \
           x.isalnum() and x or \
           len(set(string.punctuation).intersection(x)) == 1 and \
           x.count('.') == 1 and float(x) or x

## depth is wrong
## a = [1, 2, 9, 10, 17, 18, 25, 26, 33, 34, 1, 2, 3, 0, [1, 2], [8, 9], [15, 16], [22, 23]]
## M.depth(a) is 0
## sum(map(M.isiterable,a)) != len(a)

def depth(seq):
    """ for flat/nested lists find the nesting depth"""
    if not seq or not isiterable(seq): return 0
    depth = 0
    item = seq[0] # <- this is a problem: what if seq[0] is int but seq[1] is []
    previous = item
    while isiterable(item) and not isinstance(item,str):
        depth += 1
        item = item[0]
    return depth

def floatSignificantNumbers(num):
    if num == 0: return num
    sign = num > 0 and 1 or -1
    num = abs(num)
    error = 1./2**24 # 24 bits is used to store single precision float
    nlog10 = -int(N.log10(error*num))
    return sign*int(num*10**nlog10)/10.**nlog10
