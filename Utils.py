import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np

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