import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np
from dataclasses import dataclass, field
from scipy import optimize


def Kepler(E, M, e):
    return E - e * np.sin(E) - M

def volume_3d(center: tuple, dimensions: tuple) -> np.ndarray:  
    """
    Create a 3D rectangular volume centered at a specific point
    
    Parameters:
    center (tuple): (x,y,z) coordinates of center point
    dimensions (tuple): (δx, δy, δz) dimensions of the volume
    
    Returns:
    tuple: Arrays of coordinates describing the volume corners
    """
    x, y, z = center
    dx, dy, dz = dimensions
    
    # Calculate bounds
    x_min, x_max = x - dx/2, x + dx/2
    y_min, y_max = y - dy/2, y + dy/2
    z_min, z_max = z - dz/2, z + dz/2
    
    # Create corner points
    corners = np.array([
        [x_min, y_min, z_min], [x_max, y_min, z_min],
        [x_min, y_max, z_min], [x_max, y_max, z_min],
        [x_min, y_min, z_max], [x_max, y_min, z_max],
        [x_min, y_max, z_max], [x_max, y_max, z_max]
    ])
    
    return corners

# To check if a point is inside:
def is_point_in_volume(point_position: np.ndarray, point_velocity: np.ndarray, center: np.ndarray, center_v: np.ndarray, dimensions: np.ndarray, dimensions_v: np.ndarray) -> bool:
    x, y, z = point_position
    cx, cy, cz = center
    dx, dy, dz = dimensions

    vx, vy, vz = point_velocity
    c_vx, c_vy, c_vz = center_v
    d_vx, d_vy, d_vz = dimensions_v
    # print(cx, c_vx, d_vx, x, vx)
    return (abs(x - cx) <= dx/2 and 
            abs(y - cy) <= dy/2 and 
            abs(z - cz) <= dz/2 and 
            abs(vx - c_vx) <= d_vx/2 and 
            abs(vy - c_vy) <= d_vy/2 and 
            abs(vz - c_vz) <= d_vz/2)

def p_E_uniform(a_min: float, a_max: float, e_min: float, e_max: float, i_min: float, i_max: float, 
        Omega_min: float, Omega_max: float, w_min: float, w_max: float, E_min: float, E_max: float) -> float:
    p_E = 1/(a_max - a_min) * 1/(e_max - e_min) * 1/(i_max - i_min) * 1/(Omega_max - Omega_min) * 1/(w_max - w_min) * 1/(E_max - E_min)
    return p_E

#utils conversions and constants

from dataclasses import dataclass, field

@dataclass
class CanonicalUnits:
    """
    Data class to provide canonical units and compute the gravitational parameter mu in canonical units.
    """
    deg: float = field(default=np.pi / 180, init=False)
    AU_m: float = field(default=1.496e11, init=False)  # meters
    M_sun: float = field(default=1.9891e30, init=False)  # kg
    G: float = field(default=6.67430e-11, init=False)  # m^3 / (kg s^2)
    year: float = field(default=365.25 * 24 * 3600, init=False)  # seconds

    @property
    def G_cu(self) -> float:
        # Gravitational constant in canonical units
        return self.G * (1 / self.AU_m) ** 3 * self.M_sun * (self.year) ** 2

    @property
    def mu(self) -> float:
        # Gravitational parameter mu in canonical units
        return self.G_cu


from typing import Optional

@dataclass
class OrbitalElements:
    a: float
    e: float
    i: float
    Omega: float
    w: float
    E: float | None = None
    M: float | None = None

    def __post_init__(self):
        if self.E is None and self.M is None:
            raise ValueError("Either 'E' (Eccentric Anomaly) or 'M' (Mean Anomaly) must be provided.")
        if self.M is not None and self.E is None:
            self.E = optimize.newton(Kepler, 1, args=(self.M, self.e))
        if self.E is not None and self.M is None:
            # Optionally, you could allow both, or raise an error if both are provided
            pass

@dataclass
class GravitationalParameters:
    G: float | None = None
    M: float | None = None
    m: float | None = None
    mu: float | None = None
    
    def __post_init__(self):
        """Validate and calculate missing parameters after initialization"""
        # Check if we have enough information to calculate mu
        if self.mu is None:
            if self.G is not None and self.M is not None:
                self.mu = self.G * self.M
            else:
                raise ValueError("Either 'mu' must be provided, or both 'G' and 'M' must be provided")
        
        # If we have mu but not G or M, we can't calculate them uniquely
        # This is fine - the user can access mu directly
    
    @classmethod
    def from_mu(cls, mu: float):
        """Create instance with only mu parameter"""
        return cls(mu=mu)
    
    @classmethod
    def from_GM(cls, G: float, M: float, m: float | None = None):
        """Create instance with G and M parameters, optionally m"""
        return cls(G=G, M=M, m=m)
    
    @classmethod
    def from_standard(cls, m: float | None = None):
        """Create instance with standard gravitational constant and solar mass"""
        return cls(G=6.67430e-11, M=1.9891e30, m=m)

class OrbitTrasformations:

    @staticmethod
    def sqrt_e(e: float) -> float:
        return (1-e**2)**0.5

    @staticmethod
    def nu(a: float, mu: float) -> float:
        return (mu*a)**0.5
    
    @staticmethod
    def r(a: float, e: float, E: float) -> float:
        return a*(1 - e*np.cos(E))

    @staticmethod
    def h(a: float, e: float, mu: float) -> float:
        return (mu*a*(1-e**2))**0.5

    @staticmethod
    def Matrix_trasformation(i: float, w: float, Omega: float) -> np.ndarray:
        Matrix = np.zeros((3, 3))

        Matrix[0, 0] = np.cos(Omega)*np.cos(w) - np.sin(Omega)*np.cos(i)*np.sin(w)
        Matrix[0, 1] = -np.cos(Omega)*np.sin(w) - np.sin(Omega)*np.cos(i)*np.cos(w)
        Matrix[0, 2] = np.sin(Omega)*np.sin(i)

        Matrix[1, 0] = np.sin(Omega)*np.cos(w) + np.cos(Omega)*np.cos(i)*np.sin(w)
        Matrix[1, 1] = -np.sin(Omega)*np.sin(w) + np.cos(Omega)*np.cos(i)*np.cos(w)
        Matrix[1, 2] = -np.cos(Omega)*np.sin(i)

        Matrix[2, 0] = np.sin(i)*np.sin(w)
        Matrix[2, 1] = np.sin(i)*np.cos(w)
        Matrix[2, 2] = np.cos(i)

        return Matrix

    @staticmethod
    def state_vector(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        M = elements.M
        q = a*(1-e)
        mu = grav_params.mu

        if M is not None:
            state_vector = spy.conics([q, e, i, w, Omega, M]+[0, mu], 0)
            return state_vector
        else:
            raise ValueError('Mean anomaly is not provided')

        """
        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        v_r = ot.nu(a, mu)/ot.r(a, e, E)

        position = a * Matrix@np.array([np.cos(E) - e, ot.sqrt_e(e)*np.sin(E), 0])
        velocity = v_r * Matrix@np.array([-np.sin(E), ot.sqrt_e(e)*np.cos(E), 0])

        return np.concatenate((position, velocity))
        """
        

class JaccobianComponents:

    """
    @staticmethod
    def partial_a(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E
        mu = grav_params.mu

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        v_r = ot.nu(a, mu)/(2*ot.r(a, e, E))

        partial_a_xyz = Matrix@np.array([np.cos(E) - e, ot.sqrt_e(e)*np.sin(E), 0])
        partial_a_vxvyvz = v_r * Matrix@np.array([np.sin(E), ot.sqrt_e(e)*np.cos(E), 0])

        return np.concatenate((partial_a_xyz, partial_a_vxvyvz))
    """

    @staticmethod
    def partial_a(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E
        mu = grav_params.mu

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        v_r = ot.nu(a, mu)/(2*ot.r(a, e, E)*a)

        partial_a_xyz = Matrix@np.array([np.cos(E) - e, ot.sqrt_e(e)*np.sin(E), 0])
        partial_a_vxvyvz = v_r * Matrix@np.array([np.sin(E), - ot.sqrt_e(e)*np.cos(E), 0])

        return np.concatenate((partial_a_xyz, partial_a_vxvyvz))

    """ 
    @staticmethod
    def partial_e(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E
        mu = grav_params.mu

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        aE = a*np.sin(E)

        a_x = a*np.sin(E)/ot.r(a, e, E) + 1
        a_y = e/ot.sqrt_e(e) - a*ot.sqrt_e(e)*np.cos(E)/ot.r(a, e, E)
        a_z = 0

        a_vx = ot.nu(a, mu)*a*np.sin(E)/(ot.r(a, e, E)**2) * (2*np.cos(E) - a*e*np.sin(E)**2/ot.r(a, e, E))
        a_vy = -ot.nu(a, mu)*e*np.cos(E)/(ot.r(a, e, E)*ot.sqrt_e(e)) + ot.nu(a, mu)*a*ot.sqrt_e(e)/(ot.r(a, e, E)**2) * (1 - 2*np.sin(E)**2 - a*e*np.cos(E)*np.sin(E)**2/ot.r(a, e, E))
        a_vz = 0

        partial_e_xyz = aE * Matrix@np.array([a_x, a_y, a_z])
        partial_e_vxvyvz = Matrix@np.array([a_vx, a_vy, a_vz])

        return np.concatenate((partial_e_xyz, partial_e_vxvyvz))
    """
    """
    @staticmethod
    def partial_e(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E
        mu = grav_params.mu

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        partial_cE = a*np.sin(E)**2/ot.r(a, e, E)
        partial_sE = a*np.cos(E)*np.sin(E)/ot.r(a, e, E)
        v_a_r = ot.nu(a, mu)*a/ot.r(a, e, E)**2
        v_r = ot.nu(a, mu)/ot.r(a, e, E)
        partial_v_r = v_a_r * (2*np.cos(E) - a*e*np.sin(E)**2/ot.r(a, e, E)) * np.sin(E)
        partial_epsilon = -e/ot.sqrt_e(e)

        e_x = partial_cE + 1
        e_y = partial_epsilon*np.sin(E) + ot.sqrt_e(e)*partial_sE
        e_z = 0

        e_vx = -partial_v_r       
        e_vy = v_a_r * partial_epsilon * (1 - 2*np.sin(E)**2 - a*e*np.cos(E)*np.sin(E)**2/ot.r(a, e, E)) + v_r * partial_epsilon * np.cos(E)
        e_vz = 0

        partial_e_xyz = -a * Matrix@np.array([e_x, e_y, e_z])
        partial_e_vxvyvz = Matrix@np.array([e_vx, e_vy, e_vz])

        return np.concatenate((partial_e_xyz, partial_e_vxvyvz))
    """

    @staticmethod
    def partial_e(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E
        mu = grav_params.mu
        ot = OrbitTrasformations()
        r = ot.r(a, e, E)
        sinE = np.sin(E)
        #sinE = (elements.y-a*(cosE-e)*C)/(ab*eps*D)
        cosE = np.cos(E) 
        eps = ot.sqrt_e(e)
        ab = np.abs(a)
        nu = ot.nu(a, mu)
        nur = nu/r

        #dX/de
        dcosEde=-a*sinE**2/r       
        dsinEde=a*cosE*sinE/r
        dnurde=(nu*a/r**2)*(cosE-(a/r)*e*sinE**2)
        depsde=-e/eps

        drAde=a*(dcosEde-1)
        drBde=ab*(depsde*sinE+eps*dsinEde)

        dvAde=-(dnurde*sinE+nur*dsinEde)
        dvBde=(dnurde*eps*cosE+nur*depsde*cosE+nur*eps*dcosEde)

        #Trigonometric function
        cosi,sini=np.cos(i), np.sin(i)
        cosw,sinw=np.cos(w), np.sin(w)
        cosW,sinW=np.cos(Omega), np.sin(Omega)

        #Components of the rotation matrix
        A=(cosW*cosw-cosi*sinW*sinw);B=(-cosW*sinw-cosw*cosi*sinW)
        C=(cosw*sinW+sinw*cosi*cosW);D=(-sinw*sinW+cosw*cosi*cosW)
        F=sinw*sini;G=cosw*sini

        Je=np.array([
            drAde*A+drBde*B,
            drAde*C+drBde*D,
            drAde*F+drBde*G,
            dvAde*A+dvBde*B,
            dvAde*C+dvBde*D,
            dvAde*F+dvBde*G,
        ])

        return Je

    @staticmethod
    def partial_i(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        Omega = elements.Omega
        state_vector = OrbitTrasformations.state_vector(elements, grav_params)

        partial_i_xyz = np.array([state_vector[2]*np.sin(Omega), -state_vector[2]*np.cos(Omega), -state_vector[0]*np.sin(Omega) + state_vector[1]*np.cos(Omega)])
        partial_i_vxvyvz = np.array([state_vector[5]*np.sin(Omega), -state_vector[5]*np.cos(Omega), -state_vector[3]*np.sin(Omega) + state_vector[4]*np.cos(Omega)])

        return np.concatenate((partial_i_xyz, partial_i_vxvyvz))

    @staticmethod
    def partial_Omega(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        state_vector = OrbitTrasformations.state_vector(elements, grav_params)

        partial_Omega_xyz = np.array([-state_vector[1], state_vector[0], 0])
        partial_Omega_vxvyvz = np.array([-state_vector[4], state_vector[3], 0])

        return np.concatenate((partial_Omega_xyz, partial_Omega_vxvyvz))

    @staticmethod
    def partial_w(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:

        i = elements.i
        Omega = elements.Omega

        state_vector = OrbitTrasformations.state_vector(elements, grav_params)

        w_x = -state_vector[1]*np.cos(i) - state_vector[2]*np.sin(i)*np.cos(Omega)
        w_y = state_vector[0]*np.cos(i) - state_vector[2]*np.sin(i)*np.sin(Omega)
        w_z = state_vector[0]*np.sin(i)*np.cos(Omega) + state_vector[1]*np.sin(i)*np.sin(Omega)

        w_vx = -state_vector[4]*np.cos(i) - state_vector[5]*np.sin(i)*np.cos(Omega)
        w_vy = state_vector[3]*np.cos(i) - state_vector[5]*np.sin(i)*np.sin(Omega)
        w_vz = state_vector[3]*np.sin(i)*np.cos(Omega) + state_vector[4]*np.sin(i)*np.sin(Omega)

        partial_w_xyz = np.array([w_x, w_y, w_z])
        partial_w_vxvyvz = np.array([w_vx, w_vy, w_vz])

        return np.concatenate((partial_w_xyz, partial_w_vxvyvz))

    
    @staticmethod
    def partial_M(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        a = elements.a
        e = elements.e
        E = elements.E
        mu = grav_params.mu

        state_vector = OrbitTrasformations.state_vector(elements, grav_params)
        n = (mu/a**3)**0.5
        factor = -(mu*a**3)**0.5/OrbitTrasformations.r(a,e,E)**3

        partial_M_xyz = (1/n) * np.array([state_vector[3], state_vector[4], state_vector[5]])
        partial_M_vxvyvz = factor * np.array([state_vector[0], state_vector[1], state_vector[2]])

        return np.concatenate((partial_M_xyz, partial_M_vxvyvz)) 


    @staticmethod
    def Jacobian(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        """Compute the complete Jacobian matrix"""
        jacobian = np.zeros((6, 6))
        
        jacobian[:, 0] = JaccobianComponents.partial_a(elements, grav_params)
        jacobian[:, 1] = JaccobianComponents.partial_e(elements, grav_params)
        jacobian[:, 2] = JaccobianComponents.partial_i(elements, grav_params)
        jacobian[:, 3] = JaccobianComponents.partial_Omega(elements, grav_params)
        jacobian[:, 4] = JaccobianComponents.partial_w(elements, grav_params)
        jacobian[:, 5] = JaccobianComponents.partial_M(elements, grav_params)
        
        a = elements.a
        e = elements.e
        q = a/(1-e)
        
        Jq=np.eye(6)
        Jq[0,0]=1/(1-e)
        Jq[0,1]=q/(1-e)**2
        jacobianq = np.matmul(jacobian,Jq)
        
        return jacobianq, jacobian

    @staticmethod
    def Jacobian_inv(elements: OrbitalElements, grav_params: GravitationalParameters) -> np.ndarray:
        return np.linalg.inv(JaccobianComponents.Jacobian(elements, grav_params))




def calcKeplerianJacobians(mu,celements,state):
    """
    Compute the Jacobian Matrix of the transformation from classical
    orbital elements (q,e,i,w,W,M) to cartesian state vector (x,y,z,x',y',z').

    Parameters:
        mu: Gravitational parameter.
        celements: Cometary elements (q,e,i,w,W,M)

    Return:

        det JXoc, where JXoc = JXoe * Jeoc, and:

            JXoe = [dx/da,dx/de,dx/di,dx/dw,dx/dW,dx/dM,
                    dy/da,dy/de,dy/di,dy/dw,dy/dW,dy/dM,
                    dz/da,dz/de,dz/di,dz/dw,dz/dW,dz/dM,
                    dx'/da,dx'/de,dx'/di,dx'/dw,dx'/dW,dx'/dM,
                    dy'/da,dy'/de',dy'/di,dy'/dw,dy'/dW,dy'/dM,
                    dz'/da,dz'/de,dz'/di',dz'/dw,dz'/dW,dz'/dM],

                    Numpy array 6x6, units compatible with mu and a.

        and:

            Jeoc = [da/dq,da/de,...]

            Jeoc = [1/(1-e),q/(1-e)**2,0,0,0,0,
                    0      ,1         ,0,0,0,0,
                    0      ,0         ,1,0,0,0,
                    0      ,0         ,0,1,0,0,
                    0      ,0         ,0,0,1,0,
                    0      ,0         ,0,0,0,1,
                    ]

            det Jeoc = 1/(1-e)

        Jacobians are used for transform pdf:

                p(c) = p(X) det(JXoc)
    """
    q,e,i,W,w,M=celements
    a=q/(1-e)

    #Orbit signature
    if e<1:
        s=+1
    elif e>1:
        s=-1
    else:
        s=0

    #Trigonometric function
    cosi,sini=np.cos(i), np.sin(i)
    cosw,sinw=np.cos(w), np.sin(w)
    cosW,sinW=np.cos(W), np.sin(W)

    #Components of the rotation matrix
    A=(cosW*cosw-cosi*sinW*sinw);B=(-cosW*sinw-cosw*cosi*sinW)
    C=(cosw*sinW+sinw*cosi*cosW);D=(-sinw*sinW+cosw*cosi*cosW)
    F=sinw*sini;G=cosw*sini

    #Primary auxiliar variables
    ab=np.abs(a)
    n=np.sqrt(mu/ab**3)
    nu=n*a**2
    eps=np.sqrt(s*(1-e**2))

    #Get cartesian coordinates
    x,y,z,vx,vy,vz=state
    r=(x**2+y**2+z**2)**0.5
    nur=nu/r

    #Eccentric anomaly as obtained from indirect information
    #From the radial equation: r = a (1-e cos E)
    cosE=(1/e)*(1-r/a)

    #From the general equation for y
    #NOTE: This is the safest way to obtain sinE without the danger of singularities
    sinE=(y-a*(cosE-e)*C)/(ab*eps*D)

    #dX/da
    Ja=np.array([x/a,y/a,z/a,-vx/(2*a),-vy/(2*a),-vz/(2*a)])

    #dX/de
    dcosEde=-s*a*sinE**2/r
    dsinEde=a*cosE*sinE/r
    dnurde=(nu*a/r**2)*(cosE-(ab/r)*e*sinE**2)
    depsde=-s*e/eps

    drAde=a*(dcosEde-1)
    drBde=ab*(depsde*sinE+eps*dsinEde)

    dvAde=-(dnurde*sinE+nur*dsinEde)
    dvBde=(dnurde*eps*cosE+nur*depsde*cosE+nur*eps*dcosEde)

    Je=np.array([
        drAde*A+drBde*B,
        drAde*C+drBde*D,
        drAde*F+drBde*G,
        dvAde*A+dvBde*B,
        dvAde*C+dvBde*D,
        dvAde*F+dvBde*G,
    ])

    #dX/di
    Ji=np.array([z*sinW,-z*cosW,-x*sinW+y*cosW,vz*sinW,-vz*cosW,-vx*sinW+vy*cosW])

    #dX/dw
    Jw=np.array([-y*cosi-z*sini*cosW,x*cosi-z*sini*sinW,sini*(x*cosW+y*sinW),
                    -vy*cosi-vz*sini*cosW,vx*cosi-vz*sini*sinW,sini*(vx*cosW+vy*sinW)])

    #dX/dW
    JW=np.array([-y,x,0,-vy,vx,0])

    #dX/dM
    JM=np.concatenate(((ab**3/mu)**0.5*np.array([vx,vy,vz]),
                        (mu*ab**3)**0.5*np.array([-x/r**3,-y/r**3,-z/r**3])))

    #Jacobian
    JX2e=np.array([Ja,Je,Ji,JW,Jw,JM]).transpose()

    #Jacobian from classical elements (a) to cometary elements (q)
    Je2c=np.eye(6)
    Je2c[0,0]=1/(1-e)
    Je2c[0,1]=q/(1-e)**2
    JX2c=np.matmul(JX2e,Je2c)

    return JX2c, JX2e

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