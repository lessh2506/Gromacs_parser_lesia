import re, os.path, os, sys, numpy as N
import itertools, textwrap
import modysa.misc as M
import logging as L

# here we collect classes representing various index files
# which are used by analysis tools

class IndexFile:

    __patterns = dict( groups = re.compile('^\s*?\[\s+?(.+)\s+?\]$', re.M),
                       line = re.compile('^[^#\n\r]+'),
                       lists = re.compile(r"""
              ([\(\)])            # handle groups (..) 
              |([ \t]+\|[ \t]+)   # column separator |
              |([ \t]+)           # spaces/tabs (field separators)
              |([0-9]+-[0-9]+)    # ranges
              |(:[0-9]+:?[0-9]*)  # loops
              |(-?[0-9]*\.[0-9]+e?[0-9]*) # floats
              |(-?[0-9]+)         # integers
                       """, re.VERBOSE),
                       colsep = re.compile('[\t ]+0[\t ]+'))

    # NOTE: order is important in lists pattern, a rule is: first more
    # specific pattern then the more general, this is why float is before 
    # integer and column sep. before field sep.

    __tags = (('group',1), ('column_sep', 2), ('field_sep',3), ('range',4),
              ('loop',5), ('float', 6), ('integer',7))

    __str_type = '_'
    __comment_chr = '#'
    # line = re.compile('^[^#\n\r]+')

    def __init__(self, groups=None, verbose=False):
        if groups is None: self.groups = {}
        elif type(groups) is dict: self.groups = groups
        self.verbose = verbose

    def reset(self):
        self.groups = {}

    def getGroups(self, group_type):
        assert group_type in [str, int]
        if group_type is str: 
            which = lambda s: s.startswith(self.__str_type)
        else: which = lambda s: not s.startswith(self.__str_type)
        return filter(which, self.groups.keys())

    def read(self, filename):
        # strip out comments (i.e. lines starting with #) and empty lines,
        # now comments might be elsewhere in a line (everything after # is 
        # omitted)
        # data = '\n'.join([l.strip() for l in open(filename,'r') \
        #                   if re.match(self.__patterns['line'],l)])
        data = []
        for l in open(filename, 'r'):
            if not re.match(self.__patterns['line'],l): continue
            cpos = l.rfind(self.__comment_chr)
            if cpos > -1: data.append(l[:cpos].strip())
            else: data.append(l.strip())
        data = '\n'.join(data)

        items = [item for item in re.split(self.__patterns['groups'], data) \
                     if item]
        groups = {}
        for i in range(0,len(items),2):
            try:
                key, value = items[i], items[i+1]
            except IndexError:
                print i, items[i:]
                raise IndexError, 'Invalid group definition found in IndexFile'
            if groups.has_key(key): groups[key].append(value)
            else: groups[key] = [value]

        for key, value in groups.items():
            # if the name of the group starts with an underscore '_'
            # the group is supposed to contain strings (names or number ranges)
            # otherwise, we expect integer numbers
            # NOTE: that previously the group name started with an asterix '*'
            if key[0] is self.__str_type: 
                group_name = key[1:]
                # value is a list with potentially multi-line strings 
                # read from various secions of IndexFile with the same
                # group_name; this is why we *join* then *split* into 
                # individual lines
                lines = [ line.strip() \
                              for line in '\n'.join(value).split('\n') \
                              if line.strip() ]
            else: 
                group_name = key
                lines = [ tuple(map(int, line.strip().split())) \
                              for line in '\n'.join(value).split('\n') \
                              if line.strip() ]
            
            assert not groups.has_key(group_name) or \
                not self.groups.has_key('%c%s' % \
                                       (self.__str_type,group_name)), \
                'Use either "%s" or "%c%s" in the index file' % \
                tuple(group_name,self.__str_type,group_name)

            # in case we read several index files 
            if self.groups.has_key(key): self.groups[key].extend(lines)
            else: self.groups[key] = lines
        
        # # number ranges present, convert them in list of tuples
        # group_ranges = [gn for gn in self.groups.keys() if gn.startswith('*')]
        # self.range(*group_ranges)

    def write(self, filename):
        assert self.groups
        assert not os.path.lexists(filename)
        ndx = open(filename, 'w')
        for name, data in self.groups.items():
            ndx.write('\n\n\n[ %s ]\n' % name)
            depth = M.depth(data)
            assert depth < 2 # FIXME: if depth > 1: flattening needed
            if depth is 0:
                buf = '\n'.join([' '.join(map(str, data))])
                buf = textwrap.fill(buf)
            elif depth is 1:
                buf = '\n'.join([' '.join(map(str,entry)) for entry in data])
            ndx.write(buf)
            
        ndx.write('\n'*2)
        ndx.close()

    def valid(self, *keys):
        """check whether the read ndx file contains all named groups
        given by names (keys)"""
        self_keys = [k[0] is self.__str_type and k[1:] or k \
                         for k in self.groups.keys()]
        req_keys = [k[0] is self.__str_type and k[1:] or k for k in keys]
        return set(req_keys).issubset(set(self_keys))

    def flat(self, group_name, full=False):
        assert self.groups.has_key(group_name), \
            'ERROR: group %s not found' % str(group_name)
        assert M.isiterable(self.groups[group_name]), \
            'ERROR: group %s contains data which are not iterable' % \
            str(group_name)
        self.groups[group_name] = M.flatten(self.groups[group_name], full=full)

    def add(self, value, group_name, columns=None):
        assert self.groups.has_key(group_name), \
            'ERROR: group %s not found' % str(group_name)
        assert self.groups[group_name], \
            'ERROR: group %s contains no data' % str(group_name)
        # FIXME: sanity check: arrays of numbers only!
        if columns is not None: 
            orig = zip(*self.groups[group_name])
            data = [(i, orig[i]) for i in columns]
        else: 
            data = [ (None, self.groups[group_name]) ] # one column assumed
        result = {}
        for i, d in data:
            d = N.array(d)
            orig_type = d.dtype.type
            d = (d + orig_type(value)).tolist()
            if i is not None: orig[i] = d
            else: self.groups[group_name] = d
        if columns is not None: self.groups[group_name] = zip(*orig)

    def lists(self, group_name, split_columns=True):
        assert self.groups.has_key(group_name), \
            'ERROR: group %s not found' % str(group_name)
        assert group_name.startswith(self.__str_type), \
            'ERROR: group name %s does not starts with %s' % \
            (str(group_name), str(self.__str_type))
        # ---
        def tokenize(chunk):
            """ parse a string into list of tokens according to 
            a given RE pattern with groups defined for different tokens
            """
            # http://effbot.org/zone/xml-scanner.htm
            pattern = self.__patterns['lists']
            result = []
            scan = pattern.scanner(chunk).match
            while True:
                m = scan()
                if not m: break
                result.append( (m.lastindex, m.group(m.lastindex)) )
            return result
        # ---
        lines = self.groups[group_name]
        # parse lines
        GRO, COL, SEP, RAN, LOO, FLO, INT = [v for k,v in self.__tags]
        result = {0: []}
        for l in lines:
            ichunk = 0
            level = 0
            buf = {level: []}
            # loop over tokens
            if self.verbose: print '[IndexFile] parsing line: %s' % l
            for code, text in tokenize(l):
                if code == GRO:
                    # group of indices
                    if text == '(': # open a group
                        level += 1
                        buf[level] = []
                    else:           # close a group
                        level -= 1
                        assert level >= 0, \
                            'improperly nested groups; line: %s' % l
                        buf[level].append(M.flatten(buf[level+1]))

                elif code == SEP: continue

                elif code == COL:
                    # new column starts
                    assert level == 0, 'improperly nested groups'
                    if sum(map(M.isiterable, buf[level])) != \
                            len(buf[level]):
                        buf[level] = [ M.flatten(buf[level], full=False) ]
                    result[ichunk].append(buf[level])
                    buf = {level: []}
                    ichunk += 1
                    if not result.has_key(ichunk): result[ichunk] = []

                elif code == RAN:
                    # range of integer numbers
                    start, end = map(int, text.split('-'))
                    my_range = range(start, end+1)
                    buf[level].append(my_range)

                elif code == LOO:
                    # lists from loops
                    gj = map(int,text[1:].split(':'))
                    assert len(gj) == 2, \
                        'both n and offset needed (%s)' % text
                    n, offset = gj

                    group = buf[level].pop()

                    assert M.isiterable(group) and \
                        M.isflat(group) or type(group) is int, \
                        'group has to be FLAT iterable object ' \
                        'or an integer number but we have: ' + \
                        repr(group)

                    size = M.isiterable(group) and len(group) or 1
                    temp = N.resize([ group ]*n, (n, size))
                    temp += (N.arange(n)*offset)[:,N.newaxis]
                    temp = temp.tolist()

                    if level > 0: temp = M.flatten(temp)
                    buf[level].extend(temp)

                elif code == FLO:
                    # float numbers
                    buf[level].append(float(text))

                elif code == INT:
                    # integer numbers
                    buf[level].append(int(text))

            # storing parsed data
            assert level == 0, 'improperly nested groups; line: %s' % l
            if sum(map(M.isiterable, buf[level])) != len(buf[level]):
                buf[level] = [ M.flatten(buf[level], full=False) ]
            result[ichunk].append(buf[level])
        # return list of tuples each made of ncol lists
        # i.e. row by row
        # ---
        # beware of the fact that *result* is truncated so that
        # it fits the shortest column (in case of multirow columns);
        # this is a feature, not a flaw :)
        result = zip(*[M.flatten(result[i],full=False) \
                           for i in range(len(result))])
        # ok, done..
        self.groups[group_name[1:]] = result
        del(self.groups[group_name])

    def combinations(self, group_name, elements=None):
        # FIXME: all methods should work on a single group which
        # makes easy to adjust parameters just for a certain group
        assert self.groups.has_key(group_name), \
            'ERROR: group %s not found' % str(group_name)
        result = []
        for item in self.groups[group_name]:
            if len(item) > 1: 
                result.extend(itertools.product(*item))
                if elements is None: elements = len(item)
            else: 
                assert elements is not None, \
                    'ERROR: missing the number of elements in combinations'
                result.extend(itertools.combinations(item, elements))
        if result: self.groups[group_name] = sorted(list(set(result)))


if __name__ == '__main__':
    ndx = IndexFile()
    ndx.read('data/lipidA-ions.ndx')
    ndx.flat('atoms','moieties')
    ndx.combinations('pairs')
    print len(ndx.groups['pairs'])
