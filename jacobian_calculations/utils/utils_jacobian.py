import numpy as np
import spiceypy as spy
from tqdm import tqdm

def compute_uniform_elements(N: int) -> np.ndarray:
    """
    Compute uniform elements for N objects.
    """
    a_uniform = np.random.uniform(0, 2, N)
    e_uniform = np.random.uniform(0, 1, N)
    i_uniform = np.random.uniform(0, np.pi, N)
    Omega_uniform = np.random.uniform(0, 2*np.pi, N)
    w_uniform = np.random.uniform(0, 2*np.pi, N)
    M_uniform = np.random.uniform(0, 2*np.pi, N)
    q_uniform = a_uniform*(1-e_uniform)

    elements = np.column_stack((q_uniform, e_uniform, i_uniform, Omega_uniform, w_uniform, M_uniform))
    return elements


def compute_state_vector(elements: np.ndarray, mu: float) -> np.ndarray:
    N = len(elements)
    
    #los trasformo a vector de estado 
    state_vector_uniform = np.zeros((N, 6))
    for j in tqdm(range(N)):
        element = elements[j]
        et = 0
        state_vector_uniform[j] = spy.conics(list(element) + [0.0, mu], et)
    return state_vector_uniform


def count_particles_in_volume(state_vector: np.ndarray, center: list, size: list, idx: int) -> int:
    objsx = (abs(state_vector[:,0] - center[0]) <= size[0]) 
    objsy = (abs(state_vector[:,1] - center[1]) <= size[1]) 
    objsz = (abs(state_vector[:,2] - center[2]) <= size[2]) 
    objsvx = (abs(state_vector[:,3] - center[3]) <= size[3]) 
    objsvy = (abs(state_vector[:,4] - center[4]) <= size[4]) 
    objsvz = (abs(state_vector[:,5] - center[5]) <= size[5]) 

    objects = objsx * objsy * objsz * objsvx * objsvy * objsvz
    selected_objects = state_vector[objects][idx]

    return objects.sum(), selected_objects

def uniform_pdf(bounds):
    prob = 1.0
    for (low, high) in bounds:
        prob *= 1.0 / (high - low)
    return prob

def X2E(X,mu):
    elts=spy.oscelt(X,0,mu)
    E=elts[:6]
    return E


def computeNumericalJacobian(jfun,x,dx,**args):
    """
    Computes numerically the Jacobian matrix of a multivariate function.

    Parameters:
        jfun: multivariate function with the prototype "def jfun(x,**args)", function
        x: indepedent variables, numpy array (N).
        dx: step size of independent variables, numpy array (N).
        **args: argument of the function

    Return:
        y: dependent variables, y=jfun(x,**args)
        Jyx: Jacobian matrix:

            Jif= [dy_1/dx_1,dy_1/dx_2,...,dy_1/dx_N,
                dy_2/dx_1,dy_2/dx_2,...,dy_2/dx_N,
                                . . .
                dy_N/dx_1,dy_N/dx_2,...,dy_N/dx_N,]
    """
    N=len(x)
    J=np.zeros((N,N))
    y=jfun(x,**args)
    for i in range(N):
        for j in range(N):
            pre=[x[k] for k in range(j)]
            pos=[x[k] for k in range(j+1,N)]
            yi=lambda t:jfun(pre+[t]+pos,**args)[i]
            dyidxj=(yi(x[j]+dx[j])-yi(x[j]-dx[j]))/(2*dx[j])
            J[i,j]=dyidxj
    return y,J


def compute_phase_space_volume(size: list) -> float:
    Delta = (2*size[0]) * (2*size[1]) * (2*size[2]) * (2*size[3]) * (2*size[4]) * (2*size[5])
    return Delta


def compute_theoretical_count(N: int, PE: float, det: float, Delta: float) -> float:
    n_theoretical = N * PE * det * Delta
    return n_theoretical