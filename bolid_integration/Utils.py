import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np
from dataclasses import dataclass, field

def Geo2Eclip(lon, lat, alt, date=None, et=None, frame='ITRF93'):
    """
    Converts geodetic coordinates (latitude, longitude, altitude) of an impact 
    event on Earth to ecliptic J2000 coordinates.

    Parameters:
    ----------
    lon : float
        Geodetic longitude of the impact site [degrees].
    lat : float
        Geodetic latitude of the impact site [degrees].
    alt : float
        Altitude of the impact site above Earth's reference spheroid [km].
    date : str
        UTC date and time of the impact event in format 'YYYY-MM-DD HH:MM:SS'.

    Returns:
    -------
    r_earth_ecl : ndarray
        Cartesian coordinates [km] of the impact site in the Ecliptic J2000 frame.

    Notes:
    ------
    1. The function first converts geodetic coordinates to Earth-centered 
       Cartesian coordinates in the ITRF93 frame.
    2. Then, it applies a transformation matrix to convert from ITRF93 
       (Earth-fixed) to ECLIPJ2000 (inertial, ecliptic-based) frame.
    3. This transformation is necessary for orbital calculations, ensuring 
       the position is in an inertial reference frame.
    """
    deg = np.pi/180

    lon = lon*deg
    lat = lat*deg

    n, props = spy.bodvrd('399','RADII',3)
    RE_spice = props[0]  #Equatorial radius of the reference spheroid.
    RP_spice = props[2]  #Polar radius of the reference spheroid.
    f_spice = (RE_spice-RP_spice)/RE_spice # Flattening coefficient.
    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)
    #print("Equatorial and Polar Radios: ", RE_spice, RP_spice)

    if date: 
        et = spy.utc2et(date)  #Convert from UTC to ephemerides time
        #print("ET", et)
    
        r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)  #Convert geodetic coordinates to rectangular coordinates in the ITRF93 frame (rotante)
        #print("GeoRec", r_earth_fixed)
        M_ecl = spy.pxform(frame, 'ECLIPJ2000', et) 
        r_earth_ecl = spy.mxv(M_ecl, r_earth_fixed)  #from ITRF93 (rotante) frame to inertial frame ECLIPJ2000
    else: 
        print('et: ephemeris time is not provided')    
    return r_earth_ecl


def Geo2Rec(lon, lat, alt):
    """
    lon: (float) [°]
    lat: (float) [°]
    alt: (float) km
    date: (str) '2000-08-16 00:00:00'
    """
    deg = np.pi/180

    lon = lon*deg
    lat = lat*deg

    n, props = spy.bodvrd('399','RADII',3)
    RE_spice = props[0]
    RP_spice = props[2]
    #print("Equatorial and Polar Radios: ", RE_spice, RP_spice)
    f_spice = (RE_spice-RP_spice)/RE_spice
    #print(RE_spice, RP_spice)

    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)

    return r_earth_fixed

def Geo2Eclip2(lon, lat, alt, date):
    """
    lon: (float) [°]
    lat: (float) [°]
    alt: (float) km
    date: (str) '2000-08-16 00:00:00'
    """

    et = spy.utc2et(date)
    r_earth_fixed = Geo2Rec(lon, lat, alt)
    mx = spy.pxform('IAU_EARTH', 'ECLIPJ2000', et)
    r_earth_ecl = spy.mxv(mx, r_earth_fixed)

    return r_earth_ecl

def z_axis_rotation(x: float) -> np.ndarray:
    return np.array([[np.cos(x), -np.sin(x), 0],[np.sin(x), np.cos(x), 0],[0,0,1]])

def mag(x) -> float:
    """
    Compute the magnitude of a vector.
    
    Parameters:
    -----------
    x : array-like
        Input vector (can be list, tuple, numpy array, etc.)
        
    Returns:
    --------
    float
        Magnitude of the vector
    """
    return (x@x)**0.5

def get_velocity_ecliptic(vx: float, vy: float, vz: float, lon: float, lat: float, alt: float, 
                         date: str | None = None, et: float | None = None) -> np.ndarray:
    """
    Convert velocity vector from Earth-fixed to ecliptic J2000 coordinates.
    
    This function takes a velocity vector in Earth-fixed coordinates and converts it
    to ecliptic J2000 coordinates, accounting for Earth's rotation.
    
    Parameters:
    -----------
    vx, vy, vz : float
        Velocity components in Earth-fixed coordinates [km/s]
    lon, lat : float
        Geodetic longitude and latitude of the observation point [degrees]
    alt : float
        Altitude above Earth's reference spheroid [km]
    date : str, optional
        UTC date and time in format 'YYYY-MM-DD HH:MM:SS'
        Must be provided if et is not provided
    et : float, optional
        Ephemeris time (ET) in seconds past J2000
        Must be provided if date is not provided
        
    Returns:
    --------
    v_eclip : np.ndarray
        Velocity vector in ecliptic J2000 coordinates [km/s]
        
    Raises:
    -------
    ValueError
        If neither date nor et is provided, or if both are provided
        
    Notes:
    ------
    The function accounts for Earth's rotation by adding the cross product
    of Earth's angular velocity and the position vector to the input velocity.
    """
    # Input validation
    if date is not None and et is not None:
        raise ValueError("Provide either 'date' or 'et', not both")
    if date is None and et is None:
        raise ValueError("Either 'date' or 'et' must be provided")
    
    # Convert velocity to numpy array
    v = np.array([vx, vy, vz])
    
    # Get position vector in Earth-fixed coordinates
    r = Geo2Rec(lon, lat, alt)
    
    # Earth's rotation parameters
    t_sidereal = 86164.09053083288  # Sidereal day in seconds
    w_earth = 2 * np.pi / t_sidereal  # Earth's angular velocity [rad/s]
    omega = np.array([0, 0, w_earth])
    
    # Add Earth's rotation contribution to velocity
    v_E = v + spy.vcrss(omega, r)  # Velocity in Earth-fixed frame [km/s]
    
    # Convert ephemeris time if date is provided
    if date is not None:
        et = spy.utc2et(date)
    
    # Transform from Earth-fixed to ecliptic J2000 coordinates
    # At this point, et is guaranteed to be a float due to input validation
    mx = spy.pxform('ITRF93', 'ECLIPJ2000', et)  # type: ignore
    v_eclip = spy.mxv(mx, v_E)
    
    return v_eclip

def change_coord(x):
    #funcion para trasformar el formato de coordenadas terrestres que da CNEOS
    if x[-1] == 'N' or x[-1] == 'E':
        new = float(x[:-1])
    elif x[-1] == 'S' or x[-1] == 'W':
        new = -float(x[:-1])
    return new  

def mean_anomaly(e, E=False, F=False):
    if E:
        return E - e*np.sin(E)
    else: 
        E = 2 * np.arctan2(np.tan(F/2), ((1 + e)/(1 - e))**(-0.5))
        return E - e*np.sin(E)
    
def size(E, v, rho):
    return (12*E/(np.pi*rho*v**2))**(1/3)


def phi(n, alpha):
    if n == 1:
        An = 3.332
        Bn = 0.631
        Cn = 0.986
    else:
        An = 1.862
        Bn = 1.218
        Cn = 0.238

    W = np.exp(-90.56 * np.tan(alpha / 2)**2)
    phi_ns = 1 - (Cn / (0.119 + (1.1341 * np.sin(alpha)) - (0.754 * np.sin(alpha)**2)))
    phi_nl = np.exp(-An * np.tan(alpha/2)**Bn)
    return W * phi_ns + (1 - W)*phi_nl

def H_red(alpha, H):
    G = 0.15
    #print('phi_1', phi(n=1, alpha=alpha))
    #print('phi_1', phi(n=2, alpha=alpha))
    return H - 2.5*np.log10((1 - G) * phi(n=1, alpha=alpha) + G * phi(n=2, alpha=alpha))

def H_abs(p, D):
    H = -(1/0.2) * np.log10((D * p**0.5) / 1329.22) 
    return H 

def V(alpha, r, Delta, H):
    return H_red(alpha, H) + 5 * np.log10(r * Delta)

def H_red2(alpha, H):
    return H * (1 - alpha/180)

def V2(alpha, r, Delta, H):
    return H_red2(alpha, H) + 5 * np.log10(r * Delta)

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


@dataclass
class OrbitalElements:
    a: float 
    e: float 
    i: float 
    Omega: float 
    w: float 
    E: float

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
        E = elements.E
        mu = grav_params.mu

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        v_r = ot.nu(a, mu)/ot.r(a, e, E)

        position = a * Matrix@np.array([np.cos(E) - e, ot.sqrt_e(e)*np.sin(E), 0])
        velocity = v_r * Matrix@np.array([-np.sin(E), ot.sqrt_e(e)*np.cos(E), 0])

        return np.concatenate((position, velocity))


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

        partial_Omega_xyz = np.array([-state_vector[1], -state_vector[0], 0])
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
        i = elements.i
        w = elements.w
        Omega = elements.Omega
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
        
        return jacobian

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

    return JX2c

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