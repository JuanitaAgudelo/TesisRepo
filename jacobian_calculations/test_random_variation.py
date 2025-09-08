import numpy as np
from scipy.optimize import nnls
from tqdm import tqdm
import spiceypy as spy
from Utils import  surface_integral_P_X, trasformation_E_to_X, CanonicalUnits, GravitationalParameters, surface_integral_P_xyvxvy, trasformation_aewE_to_xyvxvy, trasformation_xyvxvy_to_aewE
import pandas as pd

deg = np.pi/180
AU_m = 1.496e11 #m
M_sun = 1.9891e30
G = 6.67430e-11 # m^3 / (kg s^2)
year = 365.25*24*3600 #s
mu = CanonicalUnits().mu
grav_params = GravitationalParameters(mu=mu)

# Define center and widths of the phase-space hypercube
c_x = 1
c_y = 0
#c_z = 0
v_x = 0
v_y = (mu/1)**0.5
#v_z = 0

dxyz = 0.5
dvxyz = 5000 * (1/AU_m) * year

#center = (c_x, c_y, c_z, v_x, v_y, v_z)
#widths = (dxyz, dxyz, dxyz, dvxyz, dvxyz, dvxyz)
center = (c_x, c_y, v_x, v_y)
widths = (dxyz, dxyz, dvxyz, dvxyz)

N = int(1e5)
iterations = 20
#df_count_objects = pd.DataFrame(columns=["transformation_type", "numerical_count", "theoretical_count", "objects_generated"])
df_count_objects = pd.read_csv("count_objects.csv")
list_count_objects = []

for i in tqdm(range(iterations)):
    print(f"Iteration {i}")
    print(f"Generating elements")
    a_uniform = np.random.uniform(0, 2, N)
    e_uniform = np.random.uniform(0, 1, N)
    #i_uniform = np.random.uniform(0, np.pi, N)
    #Omega_uniform = np.random.uniform(0, 2*np.pi, N)
    w_uniform = np.random.uniform(0, 2*np.pi, N)
    M_uniform = np.random.uniform(0, 2*np.pi, N)
    q_uniform = a_uniform*(1-e_uniform)

    print("Generation state vector")
    #state_vectors = np.zeros((N, 6))
    state_vectors = np.zeros((N, 4))
    for el in tqdm(range(N)):
        a = a_uniform[el]
        e = e_uniform[el]
        #i = i_uniform[el]
        #Omega = Omega_uniform[el]
        w = w_uniform[el]
        M = M_uniform[el]

        #x, y, z, vx, vy, vz = trasformation_E_to_X(a, e, i, Omega, w, M, mu)
        x, y, vx, vy = trasformation_aewE_to_xyvxvy(a, e, w, M, mu)
        #state_vectors[el] = np.array([x, y, z, vx, vy, vz])
        state_vectors[el] = np.array([x, y, vx, vy])


    # Theoretical: compute expected number of objects in the volume by integrating the distribution
        N_theoretical = surface_integral_P_xyvxvy(center, widths, n_points=8, mu=mu)
        #print(f"Theoretical (integral) number of objects in volume: {N_theoretical * N}")

        # Numerical: count number of objects in the volume from the sample
        objsx = (abs(state_vectors[:,0] - c_x) <= dxyz/2) 
        objsy = (abs(state_vectors[:,1] - c_y) <= dxyz/2) 
        objsvx = (abs(state_vectors[:,2] - v_x) <= dvxyz/2) 
        objsvy = (abs(state_vectors[:,3] - v_y) <= dvxyz/2) 

        objects = objsx & objsy & objsvx & objsvy
        N_numeric = objects.sum()

    list_count_objects.append({"transformation_type": "uniform_aewM", "numerical_count": N_numeric, "theoretical_count": N_theoretical * N, "objects_generated": N})

df = pd.concat([df_count_objects, pd.DataFrame(list_count_objects)], ignore_index=True)

df.to_csv("count_objects.csv", index=False)

