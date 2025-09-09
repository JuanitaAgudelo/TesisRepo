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


def compute_functions(i: float, w: float, Omega: float, 
                    which: Set[str]) -> Dict[str, float]:
    """Compute transformation functions A, B, C, D, F, G"""
    results = {}
    
    if 'A' in which:
        results['A'] = (np.cos(Omega)*np.cos(w) - np.sin(Omega)*np.cos(i)*np.sin(w))
    if 'B' in which:
        results['B'] = (-np.cos(Omega)*np.sin(w) - np.sin(Omega)*np.cos(i)*np.cos(w))
    if 'C' in which:
        results['C'] = (np.sin(Omega)*np.cos(w) + np.cos(Omega)*np.cos(i)*np.sin(w))
    if 'D' in which:
        results['D'] = (-np.sin(Omega)*np.sin(w) + np.cos(Omega)*np.cos(i)*np.cos(w))
    
    return results

def compute_state_vector(E: np.array, mu: float) -> np.array:
    state_vector = spy.conics(E+[0, mu], 0)
    return state_vector

def jacobian_XoE(a: float, e: float, w: float, M: float, mu: float) -> np.array:
    i = 0
    Omega = 0
    E = optimize.newton(Kepler, M, args=(M, e))

    functions = compute_functions(i, w, Omega, {'A', 'B', 'C', 'D'})
    A = functions['A']
    B = functions['B']
    C = functions['C']
    D = functions['D']
    
    ot = OrbitTrasformations()
    r = ot.r(a, e, E)
    eps = ot.sqrt_e(e)
    nu = ot.nu(a, mu)
    nur = nu/r

    q = a*(1-e)
    state_vector = compute_state_vector([q, e, i, Omega, w, M], mu)

    #partial a

    partial_a_x = state_vector[0]/a
    partial_a_y = state_vector[1]/a
    partial_a_vx = -state_vector[3]/(2*a)
    partial_a_vy = -state_vector[4]/(2*a)

    partial_a = [partial_a_x, partial_a_y, partial_a_vx, partial_a_vy]

    #partial e

    #dX/de
    dcosEde=-a*np.sin(E)**2/r       
    dsinEde=a*np.cos(E)*np.sin(E)/r
    dnurde=(nu*a/r**2)*(np.cos(E)-(a/r)*e*np.sin(E)**2)
    depsde=-e/eps

    drAde=a*(dcosEde-1)
    drBde=a*(depsde*np.sin(E)+eps*dsinEde)

    dvAde=-(dnurde*np.sin(E)+nur*dsinEde)
    dvBde=(dnurde*eps*np.cos(E)+nur*depsde*np.cos(E)+nur*eps*dcosEde)

    partial_e = np.array([
        drAde*A+drBde*B,
        drAde*C+drBde*D,
        dvAde*A+dvBde*B,
        dvAde*C+dvBde*D
    ])

    #partial w

    partial_w_x = -state_vector[1]*np.cos(i) - state_vector[2]*np.sin(i)*np.cos(Omega)
    partial_w_y = state_vector[0]*np.cos(i) - state_vector[2]*np.sin(i)*np.sin(Omega)
    partial_w_vx = -state_vector[4]*np.cos(i) - state_vector[5]*np.sin(i)*np.cos(Omega)
    partial_w_vy = state_vector[3]*np.cos(i) - state_vector[5]*np.sin(i)*np.sin(Omega)

    partial_w = [partial_w_x, partial_w_y, partial_w_vx, partial_w_vy]

    #partial M
    n = (mu/a**3)**0.5
    factor = -(mu*a**3)**0.5/r**3

    partial_M_x = (1/n) * state_vector[3]
    partial_M_y = (1/n) * state_vector[4]
    partial_M_vx = factor * state_vector[0]
    partial_M_vy = factor * state_vector[1]

    partial_M = [partial_M_x, partial_M_y, partial_M_vx, partial_M_vy]

    J = np.zeros((4,4))
    J[:,0] = partial_a
    J[:,1] = partial_e
    J[:,2] = partial_w
    J[:,3] = partial_M
    
    return J

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

def P_aewM() -> float:
    a_max = 2; a_min = 0
    e_max = 1; e_min = 0
    w_max = 2*np.pi; w_min = 0
    M_max = 2*np.pi; M_min = 0
    return 1/(a_max - a_min) * 1/(e_max - e_min) * 1/(w_max - w_min) * 1/(M_max - M_min)

def P_xyvxvy(x: float, y: float, vx: float, vy: float) -> float:
    a, e, w, M = trasformation_xyvxvy_to_aewE(x, y, vx, vy, mu)
    J = jacobian_XoE(a,e,w,M,mu)
    #det = np.linalg.det(J)
    det = 1.0/np.linalg.det(J)
    P = P_aewM() * abs(det)
    return P

def P_xyvxvy_vectorized(x: np.array, y: np.array, vx: np.array, vy: np.array, mu: float) -> np.array:
    """
    Vectorized version: x, y, vx, vy are arrays (or scalars).
    Returns array of P values.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    vx = np.asarray(vx)
    vy = np.asarray(vy)

    # Prepare output array
    shape = np.broadcast(x, y, vx, vy).shape
    P = np.empty(shape, dtype=float)

    # Flatten for iteration if needed
    x_flat = x.ravel()
    y_flat = y.ravel()
    vx_flat = vx.ravel()
    vy_flat = vy.ravel()

    for idx in range(x_flat.size):
        a, e, w, M = trasformation_xyvxvy_to_aewE(x_flat[idx], y_flat[idx], vx_flat[idx], vy_flat[idx], mu)
        J = jacobian_XoE(a,e,w,M,mu)
        #det = np.linalg.det(J)/(1-e)
        det = 1.0/np.linalg.det(J)
        P.flat[idx] = P_aewM() * abs(det)

    return P.reshape(shape)

def surface_integral_P_xyvxvy(center, widths, n_points=8, mu=1):
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

    x0, y0, vx0, vy0 = center
    dx, dy, dvx, dvy = widths

    # Get Gauss-Legendre points and weights for [-1, 1]
    pts, wts = leggauss(n_points)

    # Map points from [-1, 1] to [center-width/2, center+width/2] for each dimension
    x_pts = x0 + 0.5*dx*pts
    y_pts = y0 + 0.5*dy*pts
    vx_pts = vx0 + 0.5*dvx*pts
    vy_pts = vy0 + 0.5*dvy*pts

    # Create meshgrid of all quadrature points
    X, Y, VX, VY = np.meshgrid(x_pts, y_pts, vx_pts, vy_pts, indexing='ij')
    WX, WY, WVX, WVY = np.meshgrid(wts, wts, wts, wts, indexing='ij')

    # Flatten for vectorized evaluation
    Xf = X.ravel()
    Yf = Y.ravel()
    VXf = VX.ravel()
    VYf = VY.ravel()
    WF = (WX * WY * WVX * WVY).ravel()

    # Evaluate P at all points
    Pf = P_xyvxvy_vectorized(Xf, Yf, VXf, VYf, mu)

    # Integral is sum(P * weight) * volume factor
    integral = np.sum(Pf * WF) * (0.5*dx) * (0.5*dy) * (0.5*dvx) * (0.5*dvy)
    return integral


# Define center and widths of the phase-space hypercube
N = int(1e6)

c_x = 1
c_y = 0
v_x = 0
v_y = (mu/1)**0.5

dxy = 0.5
dvxy = 5000 * (1/AU_m) * year

center = (c_x, c_y, v_x, v_y)
widths = (dxy, dxy, dvxy, dvxy)

# Theoretical: compute expected number of objects in the volume by integrating the distribution
N_theoretical = surface_integral_P_xyvxvy(center, widths, n_points=8, mu=mu)
print(f"Theoretical (integral) number of objects in volume: {N_theoretical * N}")
# Example: Compare theoretical (integral) and numerical (count) number of objects in a phase-space volume

iterations = 18
#df_count_objects = pd.DataFrame(columns=["transformation_type", "numerical_count", "theoretical_count", "objects_generated"])
df_count_objects = pd.read_csv("count_objects.csv")
list_count_objects = []

for i in tqdm(range(iterations)):
    print(f"Iteration {i}")
    print(f"Generating elements")
    a_uniform = np.random.uniform(0, 2, N)
    e_uniform = np.random.uniform(0, 1, N)
    w_uniform = np.random.uniform(0, 2*np.pi, N)
    M_uniform = np.random.uniform(0, 2*np.pi, N)
    q_uniform = a_uniform*(1-e_uniform)

    print("Generation state vector")
    #state_vectors = np.zeros((N, 6))
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

    list_count_objects.append({"transformation_type": "uniform_aewM", "numerical_count": N_numeric, "theoretical_count": N_theoretical * N, "objects_generated": N})

df = pd.concat([df_count_objects, pd.DataFrame(list_count_objects)], ignore_index=True)

df.to_csv("count_objects.csv", index=False)

