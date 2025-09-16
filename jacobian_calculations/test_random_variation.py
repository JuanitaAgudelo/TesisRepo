import numpy as np
from scipy.optimize import nnls
from tqdm import tqdm
import spiceypy as spy
from Utils import OrbitTrasformations, Kepler, GravitationalParameters, CanonicalUnits
import pandas as pd
from typing import Set, Dict
import scipy.optimize as optimize


deg = np.pi/180
AU_m = 1.496e11 #m
M_sun = 1.9891e30
G = 6.67430e-11 # m^3 / (kg s^2)
year = 365.25*24*3600 #s
mu = CanonicalUnits().mu
grav_params = GravitationalParameters(mu=mu)


def compute_state_vector(E: np.array, mu: float) -> np.array:
    state_vector = spy.conics(E+[0, mu], 0)
    return state_vector


def trasformation_aewE_to_xyvxvy(a: float, e: float, w: float, M: float, mu: float) -> tuple[float, float, float, float]:
    Omega = 0
    i = 0
    q = a*(1-e)

    state_vec = compute_state_vector([q, e, i, Omega, w, M], mu)
    x = state_vec[0]
    y = state_vec[1]
    vx = state_vec[3]
    vy = state_vec[4]

    return x, y, vx, vy

def trasformation_xyvxvy_to_aewE(x: float, y: float, vx: float, vy: float, mu: float) -> tuple[float, float, float, float]:

    elements = spy.oscelt([x, y, 0, vx, vy, 0], et=0, mu=mu)
    q = elements[0]
    e = elements[1]
    w = elements[4]
    M = elements[5]
    a = q/(1-e)

    return a, e, w, M


# Define center and widths of the phase-space hypercube
N = int(1e7)

c_x = 1
c_y = 0
v_x = 0
v_y = (mu/1)**0.5

dxy = 0.5
dvxy = 5000 * (1/AU_m) * year

center = (c_x, c_y, v_x, v_y)
widths = (dxy, dxy, dvxy, dvxy)

iterations = 20
#df_count_objects = pd.DataFrame(columns=["transformation_type", "numerical_count", "theoretical_count", "objects_generated"])
df_count_objects = pd.read_csv("count_objects.csv")
list_count_objects = []

for i in range(iterations):
    print(f"Iteration {i}")
    a_uniform = np.random.uniform(0, 2, N)
    e_uniform = np.random.uniform(0, 1, N)
    w_uniform = np.random.uniform(0, 2*np.pi, N)
    M_uniform = np.random.uniform(0, 2*np.pi, N)
    q_uniform = a_uniform*(1-e_uniform)

    xyvxvy = np.zeros((N, 4))
    for el in tqdm(range(N)):
        a = a_uniform[el]
        e = e_uniform[el]
        w = w_uniform[el]
        M = M_uniform[el]

        #x, y, z, vx, vy, vz = trasformation_E_to_X(a, e, i, Omega, w, M, mu)
        x, y, vx, vy = trasformation_aewE_to_xyvxvy(a, e, w, M, mu)
        xyvxvy[el] = np.array([x, y, vx, vy])

    # Numerical: count number of objects in the volume from the sample
    objsx = (abs(xyvxvy[:,0] - c_x) <= dxy/2) 
    objsy = (abs(xyvxvy[:,1] - c_y) <= dxy/2) 
    objsvx = (abs(xyvxvy[:,2] - v_x) <= dvxy/2) 
    objsvy = (abs(xyvxvy[:,3] - v_y) <= dvxy/2) 

    objects = objsx & objsy & objsvx & objsvy
    N_numeric = objects.sum()

    list_count_objects.append({"transformation_type": "uniform_aewM", "numerical_count": N_numeric, "theoretical_count": 117.758252, "objects_generated": N})

df = pd.concat([df_count_objects, pd.DataFrame(list_count_objects)], ignore_index=True)

df.to_csv("count_objects.csv", index=False)

