import numpy as np
from scipy.optimize import nnls
from tqdm import tqdm
import spiceypy as spy
import sys
from pathlib import Path
import pandas as pd
from typing import Set, Dict
import scipy.optimize as optimize
sys.path.append(str(Path().resolve().parent))
from utils.Utils import (
    compute_functions, 
    compute_state_vector, 
    Kepler, OrbitTrasformations, 
    compute_jacobian_XoE,
    CanonicalUnits, 
    GravitationalParameters, 
    trasformation_E_to_X, 
    trasformation_X_to_E,
    P_E, 
    P_X_vectorized,
    )

#constants
deg = np.pi/180
AU_m = 1.496e11 #m
M_sun = 1.9891e30
G = 6.67430e-11 # m^3 / (kg s^2)
year = 365.25*24*3600 #s
mu = CanonicalUnits().mu
grav_params = GravitationalParameters(mu=mu)

def P_E_sini(i: float) -> float:
    a_max = 2; a_min = 0
    e_max = 1; e_min = 0
    
    pi = 1/2 * np.abs(np.cos(i))

    Omega_max = 2*np.pi; Omega_min = 0
    w_max = 2*np.pi; w_min = 0
    M_max = 2*np.pi; M_min = 0

    return 1/(a_max - a_min) * 1/(e_max - e_min) * pi * 1/(Omega_max - Omega_min) * 1/(w_max - w_min) * 1/(M_max - M_min)

def P_X_vectorized_sini(x: np.array, y: np.array, z: np.array, vx: np.array, vy: np.array, vz: np.array, mu: float) -> np.array:
    """
    Vectorized version: x, y, vx, vy are arrays (or scalars).
    Returns array of P values.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    vx = np.asarray(vx)
    vy = np.asarray(vy)
    vz = np.asarray(vz)
    # Prepare output array
    shape = np.broadcast(x, y, z, vx, vy, vz).shape
    P = np.empty(shape, dtype=float)

    # Flatten for iteration if needed
    x_flat = x.ravel()
    y_flat = y.ravel()
    z_flat = z.ravel()
    vx_flat = vx.ravel()
    vy_flat = vy.ravel()
    vz_flat = vz.ravel()

    for idx in range(x_flat.size):
        a, e, i, Omega, w, M = trasformation_X_to_E(x_flat[idx], y_flat[idx], z_flat[idx], vx_flat[idx], vy_flat[idx], vz_flat[idx], mu)
        J = compute_jacobian_XoE(a,e,i,Omega,w,M,mu)
        det = np.linalg.det(J) 
        inv_det = 1.0/det
        P.flat[idx] = P_E_sini(i) * abs(inv_det)

    return P.reshape(shape)

def surface_integral_P_X_sini(center, widths, n_points=8, mu=1):
    """
    Calculate the surface integral of P_xyvxvy in a hypercube centered at (x, y, vx, vy)
    with dimensions (dx, dy, dvx, dvy) using Gauss-Legendre quadrature.

    Parameters:
        center: tuple/list/array of (x, y, vx, vy) center
        widths: tuple/list/array of (dx, dy, dvx, dvy) side lengths
        n_points: number of quadrature points per dimension

    Returns:
        Integral (float)
    """
    from numpy.polynomial.legendre import leggauss

    x0, y0, z0, vx0, vy0, vz0 = center
    dx, dy, dz, dvx, dvy, dvz = widths

    # Get Gauss-Legendre points and weights for [-1, 1]
    pts, wts = leggauss(n_points)

    # Map points from [-1, 1] to [center-width/2, center+width/2] for each dimension
    x_pts = x0 + 0.5*dx*pts
    y_pts = y0 + 0.5*dy*pts
    z_pts = z0 + 0.5*dz*pts
    vx_pts = vx0 + 0.5*dvx*pts
    vy_pts = vy0 + 0.5*dvy*pts
    vz_pts = vz0 + 0.5*dvz*pts

    # Create meshgrid of all quadrature points
    X, Y, Z, VX, VY, VZ = np.meshgrid(x_pts, y_pts, z_pts, vx_pts, vy_pts, vz_pts, indexing='ij')
    WX, WY, WZ, WVX, WVY, WVZ = np.meshgrid(wts, wts, wts, wts, wts, wts, indexing='ij')

    # Flatten for vectorized evaluation
    Xf = X.ravel()
    Yf = Y.ravel()
    Zf = Z.ravel()
    VXf = VX.ravel()
    VYf = VY.ravel()
    VZf = VZ.ravel()
    WF = (WX * WY * WZ * WVX * WVY * WVZ).ravel()
    # Evaluate P at all points
    Pf = P_X_vectorized_sini(Xf, Yf, Zf, VXf, VYf, VZf, mu)

    # Integral is sum(P * weight) * volume factor
    integral = np.sum(Pf * WF) * (0.5*dx) * (0.5*dy) * (0.5*dz) * (0.5*dvx) * (0.5*dvy) * (0.5*dvz)
    #return integral, y_points
    return integral


def test_random_variation():
    #define amount of objects to generate
    N = int(1e8)
    
    #define the center and widths of the phase-space hypercube
    c_x = 1
    c_y = 0
    c_z = 0
    v_x = 0
    v_y = (mu/1)**0.5
    v_z = 0

    dxyz = 0.2
    dvxyz = 5000 * (1/AU_m) * year 

    center = (c_x, c_y, c_z, v_x, v_y, v_z)
    widths = (dxyz, dxyz, dxyz, dvxyz, dvxyz, dvxyz)

    N_theoretical = surface_integral_P_X_sini(center, widths, n_points=8, mu=mu)

    iterations = 20
    #df_count_objects = pd.DataFrame(columns=["transformation_type", "numerical_count", "theoretical_count", "objects_generated"])
    df_count_objects = pd.read_csv("count_objects2.csv")
    list_count_objects = []

    for iter in range(iterations):
        print(f"Iteration {iter}")
        a_uniform = np.random.uniform(0, 2, N)
        e_uniform = np.random.uniform(0, 1, N)
        
        #maping i as sin(i)
        u = np.random.uniform(0,1,N)
        mask = u <= 0.5
        I1 = np.arcsin(2*u[mask])
        I2 = np.pi - np.arcsin(2*(1-u[~mask]))
        I = np.concatenate([I1, I2])

        Omega_uniform = np.random.uniform(0, 2*np.pi, N)
        w_uniform = np.random.uniform(0, 2*np.pi, N)
        M_uniform = np.random.uniform(0, 2*np.pi, N)
        q_uniform = a_uniform*(1-e_uniform)

        state_vectors = np.zeros((N, 6))
        for el in tqdm(range(N)):
            a = a_uniform[el]
            e = e_uniform[el]
            i = I[el]
            Omega = Omega_uniform[el]
            w = w_uniform[el]
            M = M_uniform[el]

            x, y, z, vx, vy, vz = trasformation_E_to_X(a, e, i, Omega, w, M, mu)
            state_vectors[el] = np.array([x, y, z, vx, vy, vz])

        # Numerical: count number of objects in the volume from the sample
        objsx = (abs(state_vectors[:,0] - c_x) <= dxyz/2) 
        objsy = (abs(state_vectors[:,1] - c_y) <= dxyz/2) 
        objsz = (abs(state_vectors[:,2] - c_z) <= dxyz/2) 
        objsvx = (abs(state_vectors[:,3] - v_x) <= dvxyz/2) 
        objsvy = (abs(state_vectors[:,4] - v_y) <= dvxyz/2) 
        objsvz = (abs(state_vectors[:,5] - v_z) <= dvxyz/2) 

        objects = objsx & objsy & objsz & objsvx & objsvy & objsvz
        N_numeric = objects.sum()

        list_count_objects.append({"transformation_type": "uniform_complete_sini_2", "numerical_count": N_numeric, "theoretical_count": N_theoretical * N, "objects_generated": N, "iteration": iter})
    
    df = pd.concat([df_count_objects, pd.DataFrame(list_count_objects)])
    df.to_csv("count_objects2.csv", index=False)
    

def test_random_variation_optimize():
    print("runing test_random_variation_optimize")

    # Define center and widths of the phase-space hypercube
    N = int(1e8)
    print("N = 1e8")

    c_x = 1
    c_y = 0
    c_z = 0
    v_x = 0
    v_y = (mu/1)**0.5
    v_z = 0

    dxyz = 0.2
    dvxyz = 5000 * (1/AU_m) * year 

    center = (c_x, c_y, c_z, v_x, v_y, v_z)
    widths = (dxyz, dxyz, dxyz, dvxyz, dvxyz, dvxyz)
    
    print("Computing P_X integral")
    N_theoretical = surface_integral_P_X_sini(center, widths, n_points=8, mu=mu)

    iterations = 20
    #df_count_objects = pd.DataFrame(columns=["transformation_type", "numerical_count", "theoretical_count", "objects_generated"])
    df_count_objects = pd.read_csv("count_objects2.csv")
    list_count_objects = []

    for iter in range(iterations):
        print(f"Iteration {iter}")
        a_uniform = np.random.uniform(0, 2, N)
        e_uniform = np.random.uniform(0, 1, N)
        
        #maping i as sin(i)
        u = np.random.uniform(0,1,N)
        mask = u <= 0.5
        I1 = np.arcsin(2*u[mask])
        I2 = np.pi - np.arcsin(2*(1-u[~mask]))
        I = np.concatenate([I1, I2])

        Omega_uniform = np.random.uniform(0, 2*np.pi, N)
        w_uniform = np.random.uniform(0, 2*np.pi, N)
        M_uniform = np.random.uniform(0, 2*np.pi, N)
        
        N_numeric = 0
        for el in tqdm(range(N)):
            a = a_uniform[el]
            e = e_uniform[el]
            i = I[el]
            Omega = Omega_uniform[el]
            w = w_uniform[el]
            M = M_uniform[el]

            x, y, z, vx, vy, vz = trasformation_E_to_X(a, e, i, Omega, w, M, mu)

            # Numerical: count number of objects in the volume from the sample
            if (
                abs(x - c_x) <= dxyz/2 and 
                abs(y - c_y) <= dxyz/2 and
                abs(z - c_z) <= dxyz/2 and
                abs(vx - v_x) <= dvxyz/2 and
                abs(vy - v_y) <= dvxyz/2 and
                abs(vz - v_z) <= dvxyz/2
            ):
                N_numeric += 1

        list_count_objects.append({"transformation_type": "uniform_complete_sini_2", "numerical_count": N_numeric, "theoretical_count": N_theoretical * N, "objects_generated": N, "iteration": iter})

    df = pd.concat([df_count_objects, pd.DataFrame(list_count_objects)])
    df.to_csv("count_objects2.csv", index=False)


def __init__():
    #test_random_variation()
    test_random_variation_optimize()

if __name__ == "__main__":
    __init__()