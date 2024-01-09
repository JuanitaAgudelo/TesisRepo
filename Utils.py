import pandas as pd
from astropy.time import Time
import spiceypy as spy
import numpy as np

def Geo2Eclip(lon, lat, alt, date):
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

    et = spy.utc2et(date)
    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)
    M_itrf2ecl = spy.pxform('ITRF93', 'ECLIPJ2000', et)
    r_earth_ecl = spy.mxv(M_itrf2ecl, r_earth_fixed)

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