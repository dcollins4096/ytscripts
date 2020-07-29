
from go import *
from collections import defaultdict

class plot():
    def __init__(self,fname=".",parameters={}):
        self.fname=fname
        self.parameters=parameters

class product():
    def __init__(self, name="P", regexp="reg",myglob="glob", parameters=['core_id','frame']):
        self.myglob=myglob
        self.regexp=regexp
        self.name=name
        self.parameters=parameters
        self.plots=defaultdict(lambda: list())

    def check_glob(self):
        print("check glob")
        file_list = glob.glob(self.myglob)
        print(self.myglob)
        print(file_list)

    def get_frames(self):
        file_list = glob.glob(self.myglob)
        print(self.myglob)
        for fname in file_list:
            match = self.regexp.match(fname)
            if match is None:
                print("Error with regexp match")
                raise
            mygroups = match.groups()
            params = dict(zip(self.parameters,mygroups))
            core_id = int(params['core_id'])
            myplot = plot(fname,params)
            self.plots[core_id].append(myplot)
        if len(self.parameters) > 1:
            for p in self.parameters:
                if p != 'core_id':
                    sort_key = p
                    break
            self.plots[core_id] = sorted( self.plots[core_id],key= lambda plot : plot.parameters[sort_key])

