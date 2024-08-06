
import re

class InputFile:

    def __init__(self):
        self.data = ''
        self.fields = {}

    def load(self, filename):
        self.data = open(filename,'r').read()
        kv = re.findall('(?m)^\s*(\w+?)\s*=\s*(.+)$', self.data)
        self.fields = {}
        for k, v in kv:
            nocomment = v.find('#')
            pos = nocomment > -1 and nocomment or len(v)
            self.fields[k] = v[:pos].strip()

    def entries(self, **types):
        assert self.fields
        fields = {}
        for k,v in self.fields.items():
            # return the subset of fields with known data type
            if types.has_key(k): 
                ktype = types[k]
                assert ktype in [str, float, int, tuple] # eval not allowed
                if ktype is str: 
                    # get rid of trailing quotation marks
                    v = v.strip('\'').strip('"')
                elif ktype is tuple:
                    # simplistic but avoids eval
                    v = tuple(map(float, v[1:-1].split(',')))
                else: # float
                    try: v = ktype(v)
                    except ValueError, SyntaxError:
                        print 'Error processing entry:', k
                        raise
                fields[k] = v
        return fields

    def show(self, fields):
        print '-'*60
        for k, v in sorted(fields.items()):
            print '%20s :: %s' % (k, str(v))
        print '-'*60

    def write(self, filename):
        assert self.fields or self.data
