import pandas as pd
import numpy as np
import multimin as mm
import warnings
import pickle 
import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

df = pd.read_csv("../datos/sbdb_query_results_NEOS.csv")
df["q"] = df["a"]*(1-df["e"])
neos = df[df["q"] < 1.30]
neos

data_neas = np.array(neos[["q","e","i"]])
data_neas

scales=[1.30,1.00,180.0]
udata=np.zeros_like(data_neas)
for i in range(len(data_neas)):
    udata[i]=mm.Util.tIF(data_neas[i],scales,mm.Util.f2u)
udata 