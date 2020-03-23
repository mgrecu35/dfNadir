from   numpy import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator,LogFormatter
from matplotlib.colors import BoundaryNorm

import pickle
d=pickle.load(open('convProfs.pklz','rb'))
dist=[]
tseries=[[] for j in range(6)]
tseriesNw=[[] for j in range(6)]
datL=[]
for d1 in d:
    d1_1=d1[0]
    d1_2=d1[1]
    h1=d1_1[2]
    for r in d1_2:
        dist.append(h1-r[4])
        i0=int(h1-r[4])
        if i0>=0 and i0<6:
            tseries[i0].append([d1_1[0],r[0]])
            tseriesNw[i0].append([d1_1[0],r[2]])
            if(r[4]<2.0):
                datL.append([r[0],r[1],r[2]])
