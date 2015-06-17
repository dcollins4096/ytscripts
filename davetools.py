import matplotlib.pyplot as plt
import numpy as na
import matplotlib as mpl
import types
import glob
import os.path
import tarfile


def grep(lookfor,obj):
    if isinstance(obj, types.ListType):
        my_list = obj
    else: my_list = dir(obj)
    for i in my_list:
        if lookfor.upper() in i.upper(): print i

def stat(array,strin=''):
    print '[%0.16e,%0.16e] %s %s'%(array.min(),array.max(),array.shape,strin)

