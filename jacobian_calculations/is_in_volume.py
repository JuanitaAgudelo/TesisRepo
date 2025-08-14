import numpy as np
from scipy.optimize import nnls
from tqdm import tqdm
import spiceypy as spy
from Utils import is_point_in_volume

deg = np.pi/180
AU_m = 1.496e11 #m
M_sun = 1.9891e30
G = 6.67430e-11 # m^3 / (kg s^2)
year = 365.25*24*3600 #s
G_cu = G * (1/AU_m)**3 * (M_sun) * (year)**2  
mu = G_cu

N = int(1e10)

file_state = f"state_vector_in_volume_N_1e10.dat"
file_iterations = f"objects_in_volume.dat"

x = 1  #AU
y = 0
z = 0
center = np.array([x,y,z])

delta_x = 0.1  #AU
delta_y = 0.1  #AU
delta_z = 0.1  #AU
dimensions = np.array([delta_x, delta_y, delta_z])

x = 1  #AU 
vx = 0
vy = (mu/x)**0.5
vz = 0
center_v = np.array([vx,vy,vz])

delta_vx = 5000 * (1/AU_m) * year #AU/year
delta_vy = 5000 * (1/AU_m) * year #AU/year
delta_vz = 5000 * (1/AU_m) * year #AU/year
dimensions_v = np.array([delta_vx, delta_vy, delta_vz])

object_found = 0

state_vector_in_volume = []
for iter in tqdm(range(N)):
    #print('Generating state vector')
    
    a = np.random.uniform(0, 2, 1)[0]
    e = np.random.uniform(0, 1, 1)[0]
    i = np.random.uniform(0, 180, 1)[0]*deg
    w = np.random.uniform(0, 360, 1)[0]*deg
    Omega = np.random.uniform(0, 360, 1)[0]*deg
    M = np.random.uniform(0, 360, 1)[0]*deg
    q=a*(1-e)

    elements_spice = np.array([q, e, i, Omega, w, M, 0.0, mu])
    et = 0 
    state_vector = spy.conics(elements_spice, et)

    #print('Checking if object is in volume')
    volume_object = is_point_in_volume(state_vector[:3], state_vector[3:], 
                                    center, center_v, dimensions, dimensions_v)
    
    if volume_object:
        #print('Object is in volume')
        #print(f'Object {state_vector} in volume')
        object_found += 1
        state_vector_in_volume.append(state_vector)

state_vector_in_volume = np.array(state_vector_in_volume)
np.savetxt(file_state, state_vector_in_volume)

info_object = np.array([object_found, iter, N])
print(f'Number of objects in volume: {object_found}')
np.savetxt(file_iterations, info_object)