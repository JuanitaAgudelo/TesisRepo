import pandas as pd
import numpy as np
import multimin as mm
import warnings
import pickle 
import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

# Fit parameters
Ngauss = 10
Nvars = 3
scales=[1.30,1.00,np.pi]

deg = np.pi/180

df = pd.read_csv("../datos/sbdb_query_results_NEOS.csv")
df["q"] = df["a"]*(1-df["e"])
df["i"] = df["i"]*deg
neos = df[df["q"] < 1.30]

data_neas = np.array(neos[["q","e","i"]])
udata=np.zeros_like(data_neas)
for i in range(len(data_neas)):
    udata[i]=mm.Util.tIF(data_neas[i],scales,mm.Util.f2u)

F=mm.FitCMND(Ngauss=Ngauss,Nvars=Nvars)

t = mm.Util.elTime(0)
F.fitData(udata,verbose=False,advance=1)
t = mm.Util.elTime()
print(f"-log(L)/N = {F.solution.fun/len(udata)}")

print(" ----------------- finally fitted -----------------")
print(F.cmnd)

print(" ----------------- saving fit -----------------")
F.saveFit(f"products/fit-NEAs-qei-Ng{Ngauss}-Nv{Nvars}-rad.pkl",useprefix=False)
