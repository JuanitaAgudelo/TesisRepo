import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np
from typing import Dict, Set, Optional, Union, Tuple
from dataclasses import dataclass
#import numba

# constants
G = 6.67430e-11 # m^3 / (kg s^2)
M_sun = 1.9891e30
mu = G*M_sun

@dataclass
class OrbitalElements:
    """Container for orbital elements"""
    a: float  # semi-major axis
    e: float  # eccentricity
    i: float  # inclination
    w: float  # argument of periapsis
    Omega: float  # longitude of ascending node
    E: float  # eccentric anomaly

@dataclass
class StateVector:
    """Container for Cartesian state vector"""
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float

class OrbitalTransformations:
    """Class to handle orbital transformations and calculations"""
    
    @staticmethod
    def sqrt_e(e: float) -> float:
        """Calculate sqrt(1-e^2)"""
        return (1 - e**2)**0.5
    
    @staticmethod
    def nu(a: float) -> float:
        """Calculate sqrt(mu*a)"""
        return (mu * a)**0.5
    
    @staticmethod
    def r(a: float, e: float, E: float) -> float:
        """Calculate radial distance"""
        return a * (1 - e * np.cos(E))
    
    @staticmethod
    def h(a: float, e: float) -> float:
        """Calculate specific angular momentum"""
        return (mu * a * (1 - e**2))**0.5
    
    @staticmethod
    def compute_transformation_functions(i: float, w: float, Omega: float, 
                                       which: Set[str]) -> Dict[str, float]:
        """Compute transformation functions A, B, C, D, F, G"""
        results = {}
        
        if 'A' in which:
            results['A'] = (np.cos(w) * np.cos(Omega) - 
                          np.cos(i) * np.sin(Omega) * np.sin(w))
        if 'B' in which:
            results['B'] = (np.sin(w) * np.cos(Omega) + 
                          np.cos(i) * np.sin(Omega) * np.cos(w))
        if 'C' in which:
            results['C'] = (np.sin(Omega) * np.cos(w) + 
                          np.cos(i) * np.cos(Omega) * np.sin(w))
        if 'D' in which:
            results['D'] = (np.sin(Omega) * np.sin(w) - 
                          np.cos(i) * np.cos(Omega) * np.cos(w))
        if 'F' in which:
            results['F'] = np.sin(w) * np.sin(i)
        if 'G' in which:
            results['G'] = np.cos(w) * np.sin(i)
        
        return results
    
    @staticmethod
    def compute_derivatives(i: float, w: float, Omega: float, 
                          respect_to: str, which: Optional[Set[str]] = None) -> Dict[str, float]:
        """Compute derivatives of transformation functions"""
        if which is None:
            which = {'A', 'B', 'C', 'D', 'F', 'G'}
        
        results = {}
        
        if respect_to == 'i':
            if 'A' in which:
                results['A'] = np.sin(i) * np.sin(Omega) * np.sin(w)
            if 'B' in which:
                results['B'] = -np.sin(i) * np.sin(Omega) * np.cos(w)
            if 'C' in which:
                results['C'] = -np.sin(i) * np.cos(Omega) * np.sin(w)
            if 'D' in which:
                results['D'] = -np.sin(i) * np.cos(Omega) * np.cos(w)
            if 'F' in which:
                results['F'] = np.sin(w) * np.cos(i)
            if 'G' in which:
                results['G'] = np.cos(w) * np.cos(i)
        
        elif respect_to == 'w':
            if 'A' in which:
                results['A'] = (-np.cos(Omega) * np.sin(w) - 
                              np.cos(i) * np.cos(Omega) * np.cos(w))
            if 'B' in which:
                results['B'] = (np.sin(Omega) * np.cos(w) - 
                              np.cos(i) * np.sin(Omega) * np.sin(w))
            if 'C' in which:
                results['C'] = (-np.sin(Omega) * np.sin(w) + 
                              np.cos(i) * np.cos(Omega) * np.cos(w))
            if 'D' in which:
                results['D'] = (np.sin(Omega) * np.sin(w) + 
                              np.cos(i) * np.cos(Omega) * np.sin(w))
            if 'F' in which:
                results['F'] = np.cos(w) * np.sin(i)
            if 'G' in which:
                results['G'] = -np.sin(w) * np.sin(i)
        
        elif respect_to == 'Omega':
            if 'A' in which:
                results['A'] = (-np.sin(Omega) * np.cos(w) - 
                              np.cos(i) * np.cos(Omega) * np.sin(w))
            if 'B' in which:
                results['B'] = (-np.sin(Omega) * np.sin(w) + 
                              np.cos(i) * np.cos(Omega) * np.cos(w))
            if 'C' in which:
                results['C'] = (np.cos(Omega) * np.cos(w) - 
                              np.cos(i) * np.sin(Omega) * np.sin(w))
            if 'D' in which:
                results['D'] = (np.cos(Omega) * np.sin(w) + 
                              np.cos(i) * np.sin(Omega) * np.cos(w))
            if 'F' in which:
                results['F'] = 0
            if 'G' in which:
                results['G'] = 0
        else:
            raise ValueError(f'Invalid respect_to parameter: {respect_to}')
        
        return results
    
    @staticmethod
    def get_position_vector(elements: OrbitalElements) -> np.ndarray:
        """Get position vector from orbital elements"""
        functions = OrbitalTransformations.compute_transformation_functions(
            elements.i, elements.w, elements.Omega, {'A', 'B', 'C', 'D', 'F', 'G'}
        )
        
        x = (elements.a * (np.cos(elements.E) - elements.e) * functions['A'] - 
             elements.a * OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * functions['B'])
        y = (elements.a * (np.cos(elements.E) - elements.e) * functions['C'] - 
             elements.a * OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * functions['D'])
        z = (elements.a * (np.cos(elements.E) - elements.e) * functions['F'] + 
             elements.a * OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * functions['G'])
        
        return np.array([x, y, z])
    
    @staticmethod
    def get_velocity_vector(elements: OrbitalElements) -> np.ndarray:
        """Get velocity vector from orbital elements"""
        term = mu * elements.a / (OrbitalTransformations.h(elements.a, elements.e) * 
                                 OrbitalTransformations.r(elements.a, elements.e, elements.E))
        functions = OrbitalTransformations.compute_transformation_functions(
            elements.i, elements.w, elements.Omega, {'A', 'B', 'C', 'D', 'F', 'G'}
        )
        
        vx = (-term * (np.cos(elements.E) - elements.e) * functions['B'] - 
              term * (1 - elements.e**2)**0.5 * np.sin(elements.E) * functions['A'] - 
              (mu * elements.e / OrbitalTransformations.h(elements.a, elements.e)) * functions['B'])
        vy = (-term * (np.cos(elements.E) - elements.e) * functions['D'] - 
              term * (1 - elements.e**2)**0.5 * np.sin(elements.E) * functions['C'] - 
              (mu * elements.e / OrbitalTransformations.h(elements.a, elements.e)) * functions['D'])
        vz = (term * (np.cos(elements.E) - elements.e) * functions['G'] - 
              term * (1 - elements.e**2)**0.5 * np.sin(elements.E) * functions['F'] + 
              (mu * elements.e / OrbitalTransformations.h(elements.a, elements.e)) * functions['G'])
        
        return np.array([vx, vy, vz])
    
    @staticmethod
    def get_state_vector(elements: OrbitalElements) -> np.ndarray:
        """Get complete state vector from orbital elements"""
        position = OrbitalTransformations.get_position_vector(elements)
        velocity = OrbitalTransformations.get_velocity_vector(elements)
        return np.concatenate((position, velocity))

class JacobianCalculator:
    """Class to handle Jacobian matrix calculations"""
    
    @staticmethod
    def partial_derivative_a(elements: OrbitalElements, component: str) -> float:
        """Calculate partial derivative with respect to semi-major axis"""
        functions = OrbitalTransformations.compute_transformation_functions(
            elements.i, elements.w, elements.Omega, 
            {'A', 'B'} if component in ['x', 'vx'] else 
            {'C', 'D'} if component in ['y', 'vy'] else {'F', 'G'}
        )
        
        if component == 'x':
            return ((np.cos(elements.E) - elements.e) * functions['A'] - 
                   OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * functions['B'])
        elif component == 'y':
            return ((np.cos(elements.E) - elements.e) * functions['C'] - 
                   OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * functions['D'])
        elif component == 'z':
            return ((np.cos(elements.E) - elements.e) * functions['F'] + 
                   OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * functions['G'])
        elif component in ['vx', 'vy', 'vz']:
            factor = OrbitalTransformations.nu(elements.a) / (2 * OrbitalTransformations.r(elements.a, elements.e, elements.E))
            if component == 'vx':
                return (factor * np.sin(elements.E) * functions['A'] + 
                       factor * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * functions['B'])
            elif component == 'vy':
                return (factor * np.sin(elements.E) * functions['C'] + 
                       factor * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * functions['D'])
            else:  # vz
                return (factor * np.sin(elements.E) * functions['F'] - 
                       factor * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * functions['G'])
        else:
            raise ValueError(f"Invalid component: {component}")
    
    @staticmethod
    def partial_derivative_e(elements: OrbitalElements, component: str) -> float:
        """Calculate partial derivative with respect to eccentricity"""
        a_r = elements.a / OrbitalTransformations.r(elements.a, elements.e, elements.E)
        functions = OrbitalTransformations.compute_transformation_functions(
            elements.i, elements.w, elements.Omega, 
            {'A', 'B'} if component in ['x', 'vx'] else 
            {'C', 'D'} if component in ['y', 'vy'] else {'F', 'G'}
        )
        
        if component == 'x':
            return (-elements.a * (a_r * np.sin(elements.E)**2 + 1) * functions['A'] + 
                   elements.a * np.sin(elements.E) * (elements.e / OrbitalTransformations.sqrt_e(elements.e) - 
                   a_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E)) * functions['B'])
        elif component == 'y':
            return (-elements.a * (a_r * np.sin(elements.E)**2 + 1) * functions['C'] + 
                   elements.a * np.sin(elements.E) * (elements.e / OrbitalTransformations.sqrt_e(elements.e) - 
                   a_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E)) * functions['D'])
        elif component == 'z':
            return (-elements.a * (a_r * np.sin(elements.E)**2 + 1) * functions['F'] - 
                   elements.a * np.sin(elements.E) * (elements.e / OrbitalTransformations.sqrt_e(elements.e) - 
                   a_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E)) * functions['G'])
        elif component in ['vx', 'vy', 'vz']:
            nua_r = (OrbitalTransformations.nu(elements.a) * elements.a) / OrbitalTransformations.r(elements.a, elements.e, elements.E)**2
            nue_r = OrbitalTransformations.nu(elements.a) * elements.e / OrbitalTransformations.r(elements.a, elements.e, elements.E)
            
            if component == 'vx':
                term1 = -nua_r * np.sin(elements.E) * (2 * np.cos(elements.E) - a_r * elements.e * np.sin(elements.E)**2) * functions['A']
                term2 = -(nue_r * (1 / OrbitalTransformations.sqrt_e(elements.e)) * np.cos(elements.E) - 
                         nua_r * OrbitalTransformations.sqrt_e(elements.e) * (1 - 2 * np.sin(elements.E)**2 - 
                         a_r * elements.e * np.sin(elements.E)**2 * np.cos(elements.E))) * functions['B']
                return term1 + term2
            elif component == 'vy':
                term1 = -nua_r * np.sin(elements.E) * (2 * np.cos(elements.E) - a_r * elements.e * np.sin(elements.E)**2) * functions['C']
                term2 = -(nue_r * (1 / OrbitalTransformations.sqrt_e(elements.e)) * np.cos(elements.E) - 
                         nua_r * OrbitalTransformations.sqrt_e(elements.e) * (1 - 2 * np.sin(elements.E)**2 - 
                         a_r * elements.e * np.sin(elements.E)**2 * np.cos(elements.E))) * functions['D']
                return term1 + term2
            else:  # vz
                term1 = -nua_r * np.sin(elements.E) * (2 * np.cos(elements.E) - a_r * elements.e * np.sin(elements.E)**2) * functions['F']
                term2 = -(nue_r * (1 / OrbitalTransformations.sqrt_e(elements.e)) * np.cos(elements.E) - 
                         nua_r * OrbitalTransformations.sqrt_e(elements.e) * (1 - 2 * np.sin(elements.E)**2 - 
                         a_r * elements.e * np.sin(elements.E)**2 * np.cos(elements.E))) * functions['G']
                return term1 + term2
        else:
            raise ValueError(f"Invalid component: {component}")
    
    @staticmethod
    def partial_derivative_angle(elements: OrbitalElements, angle: str, component: str) -> float:
        """Calculate partial derivative with respect to angular elements (i, w, Omega)"""
        derivatives = OrbitalTransformations.compute_derivatives(
            elements.i, elements.w, elements.Omega, respect_to=angle,
            which={'A', 'B'} if component in ['x', 'vx'] else 
            {'C', 'D'} if component in ['y', 'vy'] else {'F', 'G'}
        )
        
        if component == 'x':
            return (elements.a * (np.cos(elements.E) - elements.e) * derivatives['A'] - 
                   elements.a * OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * derivatives['B'])
        elif component == 'y':
            return (elements.a * (np.cos(elements.E) - elements.e) * derivatives['C'] - 
                   elements.a * OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * derivatives['D'])
        elif component == 'z':
            return (elements.a * (np.cos(elements.E) - elements.e) * derivatives['F'] + 
                   elements.a * OrbitalTransformations.sqrt_e(elements.e) * np.sin(elements.E) * derivatives['G'])
        elif component in ['vx', 'vy', 'vz']:
            nu_r = OrbitalTransformations.nu(elements.a) / OrbitalTransformations.r(elements.a, elements.e, elements.E)
            if component == 'vx':
                return (-nu_r * np.sin(elements.E) * derivatives['A'] - 
                       nu_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * derivatives['B'])
            elif component == 'vy':
                return (-nu_r * np.sin(elements.E) * derivatives['C'] - 
                       nu_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * derivatives['D'])
            else:  # vz
                return (-nu_r * np.sin(elements.E) * derivatives['F'] + 
                       nu_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * derivatives['G'])
        else:
            raise ValueError(f"Invalid component: {component}")
    
    @staticmethod
    def partial_derivative_M(elements: OrbitalElements, component: str) -> float:
        """Calculate partial derivative with respect to mean anomaly"""
        a_r = elements.a**2 / OrbitalTransformations.r(elements.a, elements.e, elements.E)
        functions = OrbitalTransformations.compute_transformation_functions(
            elements.i, elements.w, elements.Omega, 
            {'A', 'B'} if component in ['x', 'vx'] else 
            {'C', 'D'} if component in ['y', 'vy'] else {'F', 'G'}
        )
        
        if component == 'x':
            return (-a_r * np.sin(elements.E) * functions['A'] - 
                   a_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * functions['B'])
        elif component == 'y':
            return (-a_r * np.sin(elements.E) * functions['C'] + 
                   a_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * functions['D'])
        elif component == 'z':
            return (-a_r * np.sin(elements.E) * functions['F'] + 
                   a_r * OrbitalTransformations.sqrt_e(elements.e) * np.cos(elements.E) * functions['G'])
        elif component in ['vx', 'vy', 'vz']:
            a_r = elements.a / (OrbitalTransformations.r(elements.a, elements.e, elements.E)**2)
            ae_r = (elements.a * elements.e) / OrbitalTransformations.r(elements.a, elements.e, elements.E)
            
            if component == 'vx':
                return (a_r * (ae_r * np.sin(elements.E)**2 + OrbitalTransformations.nu(elements.a) * np.cos(elements.E)) * functions['A'] + 
                       a_r * OrbitalTransformations.sqrt_e(elements.e) * (elements.a * elements.e * np.sin(elements.E) * np.cos(elements.E) + 
                       OrbitalTransformations.nu(elements.a) * np.sin(elements.E)) * functions['B'])
            elif component == 'vy':
                return (a_r * (ae_r * np.sin(elements.E)**2 + OrbitalTransformations.nu(elements.a) * np.cos(elements.E)) * functions['C'] + 
                       a_r * OrbitalTransformations.sqrt_e(elements.e) * (elements.a * elements.e * np.sin(elements.E) * np.cos(elements.E) + 
                       OrbitalTransformations.nu(elements.a) * np.sin(elements.E)) * functions['D'])
            else:  # vz
                return (a_r * (ae_r * np.sin(elements.E)**2 + OrbitalTransformations.nu(elements.a) * np.cos(elements.E)) * functions['F'] - 
                       a_r * OrbitalTransformations.sqrt_e(elements.e) * (elements.a * elements.e * np.sin(elements.E) * np.cos(elements.E) + 
                       OrbitalTransformations.nu(elements.a) * np.sin(elements.E)) * functions['G'])
        else:
            raise ValueError(f"Invalid component: {component}")
    
    @staticmethod
    def compute_jacobian(elements: OrbitalElements) -> np.ndarray:
        """Compute the complete Jacobian matrix"""
        jacobian = np.zeros((6, 6))
        components = ['x', 'y', 'z', 'vx', 'vy', 'vz']
        parameters = ['a', 'e', 'i', 'w', 'Omega', 'M']
        
        for i, component in enumerate(components):
            jacobian[i, 0] = JacobianCalculator.partial_derivative_a(elements, component)
            jacobian[i, 1] = JacobianCalculator.partial_derivative_e(elements, component)
            jacobian[i, 2] = JacobianCalculator.partial_derivative_angle(elements, 'i', component)
            jacobian[i, 3] = JacobianCalculator.partial_derivative_angle(elements, 'w', component)
            jacobian[i, 4] = JacobianCalculator.partial_derivative_angle(elements, 'Omega', component)
            jacobian[i, 5] = JacobianCalculator.partial_derivative_M(elements, component)
        
        return jacobian
    
    @staticmethod
    def compute_jacobian_inverse(elements: OrbitalElements) -> np.ndarray:
        """Compute the inverse of the Jacobian matrix"""
        return np.linalg.inv(JacobianCalculator.compute_jacobian(elements))

# Legacy function names for backward compatibility
def Jacobian_xE(a, e, i, w, Omega, E):
    """Legacy function - use JacobianCalculator.compute_jacobian instead"""
    elements = OrbitalElements(a=a, e=e, i=i, w=w, Omega=Omega, E=E)
    return JacobianCalculator.compute_jacobian(elements)

def Jacobian_Ex(a, e, i, w, Omega, E):
    """Legacy function - use JacobianCalculator.compute_jacobian_inverse instead"""
    elements = OrbitalElements(a=a, e=e, i=i, w=w, Omega=Omega, E=E)
    return JacobianCalculator.compute_jacobian_inverse(elements)

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
    #print("Equatorial and Polar Radios: ", RE_spice, RP_spice)

    if date: 
        et = spy.utc2et(date)  #Convert from UTC to ephemerides time
        #print("ET", et)
    
    if et is None:
        raise ValueError("Either date or et must be provided")
    
    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)  #Convert geodetic coordinates to rectangular coordinates in the ITRF93 frame (rotante)
    #print("GeoRec", r_earth_fixed)
    M_ecl = spy.pxform(frame, 'ECLIPJ2000', et) 
    r_earth_ecl = spy.mxv(M_ecl, r_earth_fixed)  #from ITRF93 (rotante) frame to inertial frame ECLIPJ2000
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

def z_axis_rotation(x):
    return np.array([[np.cos(x), -np.sin(x), 0],[np.sin(x), np.cos(x), 0],[0,0,1]])

def mag(x):
    return (x@x)**0.5

def get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=None, et=None):
    #v en km/s 
    #date en UTC
    v = np.array([vx, vy, vz]) 
    r = Geo2Rec(lon, lat, alt)

    t_sideral = 86164.09053083288 
    w_earth = 2 * np.pi / t_sideral 
    omega = np.array([0,0,w_earth]) #rad/s

    v_E = v + spy.vcrss(omega, r) #km/s
    #v_E, -v, mag(v_E), np.arccos((v@r_irtf)/(np.linalg.norm(v)*np.linalg.norm(r_irtf)))*180/np.pi

    if date:
        et = spy.utc2et(date)
    
    if et is None:
        raise ValueError("Either date or et must be provided")

    mx = spy.pxform('ITRF93', 'ECLIPJ2000', et)
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

# Note: Old duplicate functions have been removed. Use the new refactored classes instead.
# For backward compatibility, the legacy functions Jacobian_xE and Jacobian_Ex are still available.


# Note: Old duplicate functions have been removed. Use the new refactored classes instead.
# For backward compatibility, the legacy functions Jacobian_xE and Jacobian_Ex are still available.


# Example usage of the refactored code:
if __name__ == "__main__":
    # Example: Calculate orbital elements and Jacobian matrix
    
    # Define orbital elements
    elements = OrbitalElements(
        a=1.5,      # semi-major axis in AU
        e=0.1,      # eccentricity
        i=0.5,      # inclination in radians
        w=1.0,      # argument of periapsis in radians
        Omega=0.8,  # longitude of ascending node in radians
        E=0.3       # eccentric anomaly in radians
    )
    
    # Get position and velocity vectors
    position = OrbitalTransformations.get_position_vector(elements)
    velocity = OrbitalTransformations.get_velocity_vector(elements)
    state_vector = OrbitalTransformations.get_state_vector(elements)
    
    print("Position vector:", position)
    print("Velocity vector:", velocity)
    print("State vector:", state_vector)
    
    # Calculate Jacobian matrix
    jacobian = JacobianCalculator.compute_jacobian(elements)
    jacobian_inverse = JacobianCalculator.compute_jacobian_inverse(elements)
    
    print("Jacobian matrix shape:", jacobian.shape)
    print("Jacobian inverse shape:", jacobian_inverse.shape)
    
    # Verify that the inverse is correct
    identity_check = np.allclose(np.dot(jacobian, jacobian_inverse), np.eye(6))
    print("Jacobian inverse verification:", identity_check)
    
    # Legacy function usage (still works)
    legacy_jacobian = Jacobian_xE(elements.a, elements.e, elements.i, 
                                 elements.w, elements.Omega, elements.E)
    print("Legacy function works:", np.allclose(jacobian, legacy_jacobian))