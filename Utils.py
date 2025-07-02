import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np
import numba

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
    #print("Equatorial and Polar Radios: ", RE_spice, RP_spice)

    if date: 
        et = spy.utc2et(date)  #Convert from UTC to ephemerides time
        #print("ET", et)
    
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

def get_velocity_ecliptic(vx, vy, vz, lon, lat, alt, date=False, et=False):
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

def compute_function(i, w, Omega, which):
    try:
        results = {}
        if 'A' in which:
            results['A'] = np.cos(w)*np.cos(Omega) - np.cos(i)*np.sin(Omega)*np.sin(w)
        if 'B' in which:
            results['B'] = np.sin(w)*np.cos(Omega) + np.cos(i)*np.sin(Omega)*np.cos(w)
        if 'C' in which:
            results['C'] = np.sin(Omega)*np.cos(w) + np.cos(i)*np.cos(Omega)*np.sin(w)
        if 'D' in which:
            results['D'] = np.sin(Omega)*np.sin(w) - np.cos(i)*np.cos(Omega)*np.cos(w)
        if 'F' in which:
            results['F'] = np.sin(w)*np.sin(i)
        if 'G' in which:
            results['G'] = np.cos(w)*np.sin(i)

    except Exception as e:
        print('Error in compute_function: ', e)
    
    return results

def nu(a):
    return (mu*a)**0.5

def r(a, e, E):
    return a*(1 - e*np.cos(E))

def h(a, e):
    return (mu*a*(1-e**2))**0.5

def sqrt_e(e):
    return (1-e**2)**0.5

def get_position_vector(a, e, i, w, Omega, E):

    function = compute_function(i, w, Omega, which={'A', 'B', 'C', 'D', 'F', 'G'})

    x = a*(np.cos(E) - e)*function['A'] - a*sqrt_e(e)*np.sin(E)*function['B']
    y = a*(np.cos(E) - e)*function['C'] - a*sqrt_e(e)*np.sin(E)*function['D']
    z = a*(np.cos(E) - e)*function['F'] + a*sqrt_e(e)*np.sin(E)*function['G']

    return np.array([x, y, z])

def get_velocity_vector(a, e, i, w, Omega, E):

    term = mu*a/(h(a, e)*r(a, e, E))
    function = compute_function(i, w, Omega, which={'A', 'B', 'C', 'D', 'F', 'G'})
    vx = - term*(np.cos(E) - e)*function['B'] - term * (1-e**2)**0.5 * np.sin(E) * function['A'] - (mu*e/h(a, e)) * function['B']
    vy = - term*(np.cos(E) - e)*function['D'] - term * (1-e**2)**0.5 * np.sin(E) * function['C'] - (mu*e/h(a, e))*function['D']
    vz = term*(np.cos(E) - e)*function['G'] - term * (1-e**2)**0.5 * np.sin(E) * function['F'] + (mu*e/h(a, e)) * function['G']
    
    return np.array([vx, vy, vz])

def get_state_vector(a, e, i, w, Omega, E):

    position = get_position_vector(a, e, i, w, Omega, E)
    velocity = get_velocity_vector(a, e, i, w, Omega, E)

    return np.concatenate((position, velocity))

def get_derived_velocity_vector(a, e, i, w, Omega, E):

    nu = nu(a)
    r = r(a, e, E)

    function = compute_function(i, w, Omega, which={'A', 'B', 'C', 'D', 'F', 'G'})  
    vx = -(nu/r)*np.sin(E)*function['A'] - (nu/r)*(1-e**2)**0.5*np.cos(E)*function['B']
    vy = -(nu/r)*np.sin(E)*function['C'] - (nu/r)*(1-e**2)**0.5*np.cos(E)*function['D']
    vz = -(nu/r)*np.sin(E)*function['F'] + (nu/r)*(1-e**2)**0.5*np.cos(E)*function['G']
    
    return np.array([vx, vy, vz])

def compute_derivatives(i, w, Omega, respect_to=None, which=None):
    try:
        results = {}
        if respect_to == 'i':
            if which is None or 'A' in which:
                results['A'] = np.sin(i)*np.sin(Omega)*np.sin(w)
            if which is None or 'B' in which:
                results['B'] = -np.sin(i)*np.sin(Omega)*np.cos(w)
            if which is None or 'C' in which:
                results['C'] = -np.sin(i)*np.cos(Omega)*np.sin(w)
            if which is None or 'D' in which:
                results['D'] = -np.sin(i)*np.cos(Omega)*np.cos(w)
            if which is None or 'F' in which:
                results['F'] = np.sin(w)*np.cos(i)
            if which is None or 'G' in which:
                results['G'] = np.cos(w)*np.cos(i)
        
        elif respect_to == 'w':
            if which is None or 'A' in which:
                results['A'] = -np.cos(Omega)*np.sin(w) - np.cos(i)*np.cos(Omega)*np.cos(w)
            if which is None or 'B' in which:
                results['B'] = np.sin(Omega)*np.cos(w) - np.cos(i)*np.sin(Omega)*np.sin(w)
            if which is None or 'C' in which:
                results['C'] = -np.sin(Omega)*np.sin(w) + np.cos(i)*np.cos(Omega)*np.cos(w)
            if which is None or 'D' in which:
                results['D'] = np.sin(Omega)*np.sin(w) + np.cos(i)*np.cos(Omega)*np.sin(w)
            if which is None or 'F' in which:
                results['F'] = np.cos(w)*np.sin(i)
            if which is None or 'G' in which:
                results['G'] = -np.sin(w)*np.sin(i)
        
        elif respect_to == 'Omega':
            if which is None or 'A' in which:
                results['A'] = -np.sin(Omega)*np.cos(w) - np.cos(i)*np.cos(Omega)*np.sin(w)
            if which is None or 'B' in which:
                results['B'] = -np.sin(Omega)*np.sin(w) + np.cos(i)*np.cos(Omega)*np.cos(w)
            if which is None or 'C' in which:
                results['C'] = np.cos(Omega)*np.cos(w) - np.cos(i)*np.sin(Omega)*np.sin(w)
            if which is None or 'D' in which:
                results['D'] = np.cos(Omega)*np.sin(w) + np.cos(i)*np.sin(Omega)*np.cos(w)
            if which is None or 'F' in which:
                results['F'] = 0
            if which is None or 'G' in which:
                results['G'] = 0
        else:
            print('Error in compute_derivatives: respect_to is not valid')

    except Exception as e:
        print('Error in compute_derivatives: ', e)
    
    return results

#----------------Derivates respecto to a----------------

def partial_aX(e, i, w, Omega, E):
    function = compute_function(i, w, Omega, which={'A', 'B'}) 
    return (np.cos(E) - e)*function['A'] - sqrt_e(e)*np.sin(E)*function['B']

def partial_aY(e, i, w, Omega, E):
    function = compute_function(i, w, Omega, which={'C', 'D'})
    return (np.cos(E) - e)*function['C'] - sqrt_e(e)*np.sin(E)*function['D']

def partial_aZ(e, i, w, Omega, E):
    function = compute_function(i, w, Omega, which={'F', 'G'}) 
    return (np.cos(E) - e)*function['F'] + sqrt_e(e)*np.sin(E)*function['G']

def partial_aVx(a, e, i, w, Omega, E):
    function = compute_function(i, w, Omega, which={'A', 'B'})
    return nu(a)/(2*r(a, e, E))*np.sin(E)*function['A'] + nu(a)/(2*r(a, e, E))*sqrt_e(e)*np.cos(E)*function['B']

def partial_aVy(a, e, i, w, Omega, E):
    function = compute_function(i, w, Omega, which={'C', 'D'})
    return nu(a)/(2*r(a, e, E))*np.sin(E)*function['C'] + nu(a)/(2*r(a, e, E))*sqrt_e(e)*np.cos(E)*function['D']

def partial_aVz(a, e, i, w, Omega, E):
    function = compute_function(i, w, Omega, which={'F', 'G'})
    return nu(a)/(2*r(a, e, E))*np.sin(E)*function['F'] - nu(a)/(2*r(a, e, E))*sqrt_e(e)*np.cos(E)*function['G']


#----------------Derivates respecto e----------------

def partial_eX(a, e, i, w, Omega, E):
    a_r = a/r(a, e, E)
    function = compute_function(i, w, Omega, which={'A', 'B'})
    return -a*(a_r*np.sin(E)**2 + 1)*function['A'] + a*np.sin(E)*(e/sqrt_e(e) - a_r*sqrt_e(e)*np.cos(E))*function['B']

def partial_eY(a, e, i, w, Omega, E):
    a_r = a/r(a, e, E)
    function = compute_function(i, w, Omega, which={'C', 'D'})
    return -a*(a_r*np.sin(E)**2 + 1)*function['C'] + a*np.sin(E)*(e/sqrt_e(e) - a_r*sqrt_e(e)*np.cos(E))*function['D']

def partial_eZ(a, e, i, w, Omega, E):
    a_r = a/r(a, e, E)
    function = compute_function(i, w, Omega, which={'F', 'G'})
    return -a*(a_r*np.sin(E)**2 + 1)*function['F'] - a*np.sin(E)*(e/sqrt_e(e) - a_r*sqrt_e(e)*np.cos(E))*function['G']

def partial_eVx(a, e, i, w, Omega, E):
    a_r = a/r(a, e, E)
    nua_r = (nu(a)*a)/r(a, e, E)**2
    nue_r = nu(a)*e/r(a, e, E)

    function = compute_function(i, w, Omega, which={'A', 'B'})

    term1 = -nua_r*np.sin(E)*(2*np.cos(E) - a_r*e*np.sin(E)**2)*function['A'] 
    term2 = -(nue_r*(1/sqrt_e(e))*np.cos(E) - nua_r*sqrt_e(e)*(1 - 2*np.sin(E)**2 - a_r*e*np.sin(E)**2*np.cos(E)))*function['B']

    return term1 + term2

def partial_eVy(a, e, i, w, Omega, E):
    a_r = a/r(a, e, E)
    nua_r = (nu(a)*a)/r(a, e, E)**2
    nue_r = nu(a)*e/r(a, e, E)

    function = compute_function(i, w, Omega, which={'C', 'D'})

    term1 = - nua_r*np.sin(E)*(2*np.cos(E) - a_r*e*np.sin(E)**2)*function['C'] 
    term2 = - (nue_r*(1/sqrt_e(e))*np.cos(E) - nua_r*sqrt_e(e)*(1 - 2*np.sin(E)**2 - a_r*e*np.sin(E)**2*np.cos(E)))*function['D']

    return term1 + term2

def partial_eVz(a, e, i, w, Omega, E):
    a_r = a/r(a, e, E)
    nua_r = (nu(a)*a)/r(a, e, E)**2
    nue_r = nu(a)*e/r(a, e, E)

    function = compute_function(i, w, Omega, which={'F', 'G'})

    term1 = - nua_r*np.sin(E)*(2*np.cos(E) - a_r*e*np.sin(E)**2)*function['F'] 
    term2 = - (nue_r*(1/sqrt_e(e))*np.cos(E) - nua_r*sqrt_e(e)*(1 - 2*np.sin(E)**2 - a_r*e*np.sin(E)**2*np.cos(E)))*function['G']

    return term1 + term2


#----------------Derivates respecto i----------------

def partial_iX(a, e, i, w, Omega, E):
    partial_i = compute_derivatives(i, w, Omega, respect_to='i', which={'A', 'B'})
    return a*(np.cos(E) - e)*partial_i['A'] - a*sqrt_e(e)*np.sin(E)*partial_i['B']

def partial_iY(a, e, i, w, Omega, E):
    partial_i = compute_derivatives(i, w, Omega, respect_to='i', which={'C', 'D'})
    return a*(np.cos(E) - e)*partial_i['C'] - a*sqrt_e(e)*np.sin(E)*partial_i['D']

def partial_iZ(a, e, i, w, Omega, E):
    partial_i = compute_derivatives(i, w, Omega, respect_to='i', which={'F', 'G'})
    return a*(np.cos(E) - e)*partial_i['F'] + a*sqrt_e(e)*np.sin(E)*partial_i['G']

def partial_iVx(a, e, i, w, Omega, E):
    partial_i = compute_derivatives(i, w, Omega, respect_to='i', which={'A', 'B'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_i['A'] - nu_r*sqrt_e(e)*np.cos(E)*partial_i['B']

def partial_iVy(a, e, i, w, Omega, E):
    partial_i = compute_derivatives(i, w, Omega, respect_to='i', which={'C', 'D'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_i['C'] - nu_r*sqrt_e(e)*np.cos(E)*partial_i['D']

def partial_iVz(a, e, i, w, Omega, E):
    partial_i = compute_derivatives(i, w, Omega, respect_to='i', which={'F', 'G'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_i['F'] + nu_r*sqrt_e(e)*np.cos(E)*partial_i['G']


#----------------Derivates respecto w----------------

def partial_wX(a, e, i, w, Omega, E):
    partial_w = compute_derivatives(i, w, Omega, respect_to='w', which={'A', 'B'})
    return a*(np.cos(E) - e)*partial_w['A'] - a*sqrt_e(e)*np.sin(E)*partial_w['B']

def partial_wY(a, e, i, w, Omega, E):
    partial_w = compute_derivatives(i, w, Omega, respect_to='w', which={'C', 'D'})
    return a*(np.cos(E) - e)*partial_w['C'] - a*sqrt_e(e)*np.sin(E)*partial_w['D']

def partial_wZ(a, e, i, w, Omega, E):
    partial_w = compute_derivatives(i, w, Omega, respect_to='w', which={'F', 'G'})
    return a*(np.cos(E) - e)*partial_w['F'] + a*sqrt_e(e)*np.sin(E)*partial_w['G']

def partial_wVx(a, e, i, w, Omega, E):
    partial_w = compute_derivatives(i, w, Omega, respect_to='w', which={'A', 'B'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_w['A'] - nu_r*sqrt_e(e)*np.cos(E)*partial_w['B']

def partial_wVy(a, e, i, w, Omega, E):
    partial_w = compute_derivatives(i, w, Omega, respect_to='w', which={'C', 'D'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_w['C'] - nu_r*sqrt_e(e)*np.cos(E)*partial_w['D']

def partial_wVz(a, e, i, w, Omega, E):
    partial_w = compute_derivatives(i, w, Omega, respect_to='w', which={'F', 'G'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_w['F'] + nu_r*sqrt_e(e)*np.cos(E)*partial_w['G']


#----------------Derivates respecto Omega----------------

def partial_OmegaX(a, e, i, w, Omega, E):
    partial_Omega = compute_derivatives(i, w, Omega, respect_to='Omega', which={'A', 'B'})
    return a*(np.cos(E) - e)*partial_Omega['A'] - a*sqrt_e(e)*np.sin(E)*partial_Omega['B']

def partial_OmegaY(a, e, i, w, Omega, E):
    partial_Omega = compute_derivatives(i, w, Omega, respect_to='Omega', which={'C', 'D'})
    return a*(np.cos(E) - e)*partial_Omega['C'] - a*sqrt_e(e)*np.sin(E)*partial_Omega['D']

def partial_OmegaZ(a, e, i, w, Omega, E):
    partial_Omega = compute_derivatives(i, w, Omega, respect_to='Omega', which={'F', 'G'})
    return a*(np.cos(E) - e)*partial_Omega['F'] + a*sqrt_e(e)*np.sin(E)*partial_Omega['G']

def partial_OmegaVx(a, e, i, w, Omega, E):
    partial_Omega = compute_derivatives(i, w, Omega, respect_to='Omega', which={'A', 'B'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_Omega['A'] - nu_r*sqrt_e(e)*np.cos(E)*partial_Omega['B']

def partial_OmegaVy(a, e, i, w, Omega, E):
    partial_Omega = compute_derivatives(i, w, Omega, respect_to='Omega', which={'C', 'D'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_Omega['C'] - nu_r*sqrt_e(e)*np.cos(E)*partial_Omega['D']

def partial_OmegaVz(a, e, i, w, Omega, E):
    partial_Omega = compute_derivatives(i, w, Omega, respect_to='Omega', which={'F', 'G'})
    nu_r = nu(a)/r(a, e, E)
    return -nu_r*np.sin(E)*partial_Omega['F'] + nu_r*sqrt_e(e)*np.cos(E)*partial_Omega['G']


#----------------Derivates respecto Mean Anomaly----------------

def partial_MX(a, e, i, w, Omega, E):
    a_r = a**2/r(a, e, E)
    function = compute_function(i, w, Omega, which={'A', 'B'}) 
    return -a_r*np.sin(E)*function['A'] - a_r*sqrt_e(e)*np.cos(E)*function['B']

def partial_MY(a, e, i, w, Omega, E):
    a_r = a**2/r(a, e, E)
    function = compute_function(i, w, Omega, which={'C', 'D'}) 
    return -a_r*np.sin(E)*function['C'] + a_r*sqrt_e(e)*np.cos(E)*function['D']

def partial_MZ(a, e, i, w, Omega, E):
    a_r = a**2/r(a, e, E)
    function = compute_function(i, w, Omega, which={'F', 'G'}) 
    return -a_r*np.sin(E)*function['F'] + a_r*sqrt_e(e)*np.cos(E)*function['G']

def partial_MVx(a, e, i, w, Omega, E):
    a_r = a/(r(a, e, E)**2)
    ae_r = (a*e)/r(a, e, E)
    function = compute_function(i, w, Omega, which={'A', 'B'})
    return a_r*(ae_r*np.sin(E)**2 + nu(a)*np.cos(E))*function['A'] + a_r*sqrt_e(e)*(a*e*np.sin(E)*np.cos(E) + nu(a)*np.sin(E))*function['B']

def partial_MVy(a, e, i, w, Omega, E):
    a_r = a/(r(a, e, E)**2)
    ae_r = (a*e)/r(a, e, E)
    function = compute_function(i, w, Omega, which={'C', 'D'})
    return a_r*(ae_r*np.sin(E)**2 + nu(a)*np.cos(E))*function['C'] + a_r*sqrt_e(e)*(a*e*np.sin(E)*np.cos(E) + nu(a)*np.sin(E))*function['D']

def partial_MVz(a, e, i, w, Omega, E):
    a_r = a/(r(a, e, E)**2)
    ae_r = (a*e)/r(a, e, E)
    function = compute_function(i, w, Omega, which={'F', 'G'})
    return a_r*(ae_r*np.sin(E)**2 + nu(a)*np.cos(E))*function['F'] - a_r*sqrt_e(e)*(a*e*np.sin(E)*np.cos(E) + nu(a)*np.sin(E))*function['G']


#----------------Jaccobian----------------

def Jacobian_xE(a, e, i, w, Omega, E):
    """
    Constructs the Jacobian matrix for the transformation between orbital elements
    and Cartesian state vector.
    
    Parameters:
    a, e, i, w, Omega: Orbital elements (semi-major axis, eccentricity,
                       inclination, argument of periapsis, longitude of ascending node)
    E: Eccentric anomaly
    
    Returns:
    jacobian: 6×6 numpy array representing the Jacobian matrix where
              rows = [x, y, z, vx, vy, vz]
              columns = [a, e, i, w, Omega, M]
    """
    jacobian = np.zeros((6, 6))
    
    # First row: Partial derivatives of x with respect to orbital elements
    jacobian[0, 0] = partial_aX(e, i, w, Omega, E)
    jacobian[0, 1] = partial_eX(a, e, i, w, Omega, E)
    jacobian[0, 2] = partial_iX(a, e, i, w, Omega, E)
    jacobian[0, 3] = partial_OmegaX(a, e, i, w, Omega, E)
    jacobian[0, 4] = partial_wX(a, e, i, w, Omega, E)
    jacobian[0, 5] = partial_MX(a, e, i, w, Omega, E)
    
    # Second row: Partial derivatives of y with respect to orbital elements
    jacobian[1, 0] = partial_aY(e, i, w, Omega, E)
    jacobian[1, 1] = partial_eY(a, e, i, w, Omega, E)
    jacobian[1, 2] = partial_iY(a, e, i, w, Omega, E)
    jacobian[1, 3] = partial_OmegaY(a, e, i, w, Omega, E)
    jacobian[1, 4] = partial_wY(a, e, i, w, Omega, E)
    jacobian[1, 5] = partial_MY(a, e, i, w, Omega, E)
    
    # Third row: Partial derivatives of z with respect to orbital elements
    jacobian[2, 0] = partial_aZ(e, i, w, Omega, E)
    jacobian[2, 1] = partial_eZ(a, e, i, w, Omega, E)
    jacobian[2, 2] = partial_iZ(a, e, i, w, Omega, E)
    jacobian[2, 3] = partial_OmegaZ(a, e, i, w, Omega, E)
    jacobian[2, 4] = partial_wZ(a, e, i, w, Omega, E)
    jacobian[2, 5] = partial_MZ(a, e, i, w, Omega, E)
    
    # Fourth row: Partial derivatives of vx with respect to orbital elements
    jacobian[3, 0] = partial_aVx(a, e, i, w, Omega, E)
    jacobian[3, 1] = partial_eVx(a, e, i, w, Omega, E)
    jacobian[3, 2] = partial_iVx(a, e, i, w, Omega, E)
    jacobian[3, 3] = partial_OmegaVx(a, e, i, w, Omega, E)
    jacobian[3, 4] = partial_wVx(a, e, i, w, Omega, E)
    jacobian[3, 5] = partial_MVx(a, e, i, w, Omega, E)
    
    # Fifth row: Partial derivatives of vy with respect to orbital elements
    jacobian[4, 0] = partial_aVy(a, e, i, w, Omega, E)
    jacobian[4, 1] = partial_eVy(a, e, i, w, Omega, E)
    jacobian[4, 2] = partial_iVy(a, e, i, w, Omega, E)
    jacobian[4, 3] = partial_OmegaVy(a, e, i, w, Omega, E)
    jacobian[4, 4] = partial_wVy(a, e, i, w, Omega, E)
    jacobian[4, 5] = partial_MVy(a, e, i, w, Omega, E)
    
    # Sixth row: Partial derivatives of vz with respect to orbital elements
    jacobian[5, 0] = partial_aVz(a, e, i, w, Omega, E)
    jacobian[5, 1] = partial_eVz(a, e, i, w, Omega, E)
    jacobian[5, 2] = partial_iVz(a, e, i, w, Omega, E)
    jacobian[5, 3] = partial_OmegaVz(a, e, i, w, Omega, E)
    jacobian[5, 4] = partial_wVz(a, e, i, w, Omega, E)
    jacobian[5, 5] = partial_MVz(a, e, i, w, Omega, E)
    
    return jacobian


def Jacobian_Ex(a, e, i, w, Omega, E):
    return np.linalg.inv(Jacobian_xE(a, e, i, w, Omega, E))