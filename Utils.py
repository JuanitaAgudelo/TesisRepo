import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np

def Geo2Eclip(lon, lat, alt, date, frame='ITRF93'):
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

    et = spy.utc2et(date)  #Convert from UTC to ephemerides time
    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)  #Convert geodetic coordinates to rectangular coordinates in the ITRF93 frame (rotante)
    M_itrf2ecl = spy.pxform(frame, 'ECLIPJ2000', et) 
    r_earth_ecl = spy.mxv(M_itrf2ecl, r_earth_fixed)  #from ITRF93 (rotante) frame to inertial frame ECLIPJ2000

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
    f_spice = (RE_spice-RP_spice)/RE_spice

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