import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np
from dataclasses import dataclass

# constants
G = 6.67430e-11 # m^3 / (kg s^2)
M_sun = 1.9891e30
mu = G*M_sun

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

@dataclass
class OrbitalElements:
    a: float 
    e: float 
    i: float 
    Omega: float 
    w: float 
    E: float

class OrbitTrasformations:

    @staticmethod
    def sqrt_e(e: float) -> float:
        return (1-e**2)**0.5

    @staticmethod
    def nu(a: float) -> float:
        return (mu*a)**0.5
    
    @staticmethod
    def r(a: float, e: float, E: float) -> float:
        return a*(1 - e*np.cos(E))

    @staticmethod
    def h(a: float, e: float) -> float:
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
    def state_vector(elements: OrbitalElements) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        v_r = ot.nu(a)/ot.r(a, e, E)

        position = a * Matrix@np.array([np.cos(E) - e, ot.sqrt_e(e)*np.sin(E), 0])
        velocity = v_r * Matrix@np.array([-np.sin(E), ot.sqrt_e(e)*np.cos(E), 0])

        return np.concatenate((position, velocity))


class JaccobianComponents:

    @staticmethod
    def partial_a(elements: OrbitalElements) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        v_ra = ot.nu(a)/(2*ot.r(a, e, E)*a)

        partial_a_xyz = Matrix@np.array([np.cos(E) - e, ot.sqrt_e(e)*np.sin(E), 0])
        partial_a_vxvyvz = v_ra * Matrix@np.array([np.sin(E), np.cos(E), 0])

        return np.concatenate((partial_a_xyz, partial_a_vxvyvz))

    @staticmethod
    def partial_e(elements: OrbitalElements) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E

        Matrix = OrbitTrasformations.Matrix_trasformation(i, w, Omega)
        ot = OrbitTrasformations()
        aE = a*np.sin(E)

        a_x = a*np.sin(E)/ot.r(a, e, E) + 1
        a_y = e/ot.sqrt_e(e) - a*ot.sqrt_e(e)*np.cos(E)/ot.r(a, e, E)
        a_z = 0

        a_vx = ot.nu(a)*a*np.sin(E)/(ot.r(a, e, E)**2) * (2*np.cos(E) - a*e*np.sin(E)**2/ot.r(a, e, E))
        a_vy = -ot.nu(a)*e*np.cos(E)/(ot.r(a, e, E)*ot.sqrt_e(e)) + ot.nu(a)*a*ot.sqrt_e(e)/(ot.r(a, e, E)**2) * (1 - 2*np.sin(E)**2 - a*e*np.cos(E)*np.sin(E)**2/ot.r(a, e, E))
        a_vz = 0

        partial_e_xyz = aE * Matrix@np.array([a_x, a_y, a_z])
        partial_e_vxvyvz = Matrix@np.array([a_vx, a_vy, a_vz])

        return np.concatenate((partial_e_xyz, partial_e_vxvyvz))

    @staticmethod
    def partial_i(elements: OrbitalElements) -> np.ndarray:
        Omega = elements.Omega

        state_vector = OrbitTrasformations.state_vector(elements)

        partial_i_xyz = np.array([state_vector[2]*np.sin(Omega), -state_vector[2]*np.cos(Omega), -state_vector[0]*np.sin(Omega) + state_vector[1]*np.cos(Omega)])
        partial_i_vxvyvz = np.array([state_vector[5]*np.sin(Omega), -state_vector[5]*np.cos(Omega), -state_vector[4]*np.sin(Omega) + state_vector[3]*np.cos(Omega)])

        return np.concatenate((partial_i_xyz, partial_i_vxvyvz))

    @staticmethod
    def partial_Omega(elements: OrbitalElements) -> np.ndarray:
        Omega = elements.Omega

        state_vector = OrbitTrasformations.state_vector(elements)

        partial_Omega_xyz = np.array([-state_vector[1]*np.sin(Omega), -state_vector[0]*np.cos(Omega), 0])
        partial_Omega_vxvyvz = np.array([-state_vector[4]*np.sin(Omega), state_vector[3]*np.cos(Omega), 0])

        return np.concatenate((partial_Omega_xyz, partial_Omega_vxvyvz))

    @staticmethod
    def partial_w(elements: OrbitalElements) -> np.ndarray:
        i = elements.i
        Omega = elements.Omega

        state_vector = OrbitTrasformations.state_vector(elements)

        w_x = -state_vector[1]*np.cos(i) - state_vector[2]*np.sin(i)*np.cos(Omega)*np.sin(Omega)
        w_y = state_vector[0]*np.sin(i) - state_vector[2]*np.sin(i)*np.sin(Omega)
        w_z = state_vector[0]*np.sin(i)*np.cos(Omega) + state_vector[1]*np.sin(i)*np.sin(Omega)

        w_vx = -state_vector[4]*np.cos(i) - state_vector[5]*np.sin(i)*np.cos(Omega)
        w_vy = state_vector[3]*np.cos(i) - state_vector[5]*np.sin(i)*np.sin(Omega)
        w_vz = state_vector[3]*np.sin(i)*np.cos(Omega) + state_vector[4]*np.sin(i)*np.sin(Omega)

        partial_w_xyz = np.array([w_x, w_y, w_z])
        partial_w_vxvyvz = np.array([w_vx, w_vy, w_vz])

        return np.concatenate((partial_w_xyz, partial_w_vxvyvz))

    
    @staticmethod
    def partial_M(elements: OrbitalElements) -> np.ndarray:
        a = elements.a
        e = elements.e
        i = elements.i
        w = elements.w
        Omega = elements.Omega
        E = elements.E

        state_vector = OrbitTrasformations.state_vector(elements)
        n = (mu/a**3)**0.5
        factor = -(mu*a**3)**0.5/OrbitTrasformations.r(a,e,E)**3

        partial_M_xyz = (1/n) * np.array([state_vector[3]*np.sin(Omega), state_vector[4]*np.cos(Omega), state_vector[5]])
        partial_M_vxvyvz = factor * np.array([state_vector[0], state_vector[1], state_vector[2]])

        return np.concatenate((partial_M_xyz, partial_M_vxvyvz)) 


    @staticmethod
    def Jacobian(elements: OrbitalElements) -> np.ndarray:
        """Compute the complete Jacobian matrix"""
        jacobian = np.zeros((6, 6))
        
        jacobian[:, 0] = JaccobianComponents.partial_a(elements)
        jacobian[:, 1] = JaccobianComponents.partial_e(elements)
        jacobian[:, 2] = JaccobianComponents.partial_i(elements)
        jacobian[:, 3] = JaccobianComponents.partial_w(elements)
        jacobian[:, 4] = JaccobianComponents.partial_Omega(elements)
        jacobian[:, 5] = JaccobianComponents.partial_M(elements)
        
        return jacobian

    @staticmethod
    def Jacobian_inv(elements: OrbitalElements) -> np.ndarray:
        return np.linalg.inv(JaccobianComponents.Jacobian(elements))

