#!/usr/bin/python

# Some tools for working with navigation data (used for SASAR)
# code: J.M. Horrell - UCT RRSG 1999-07-21

from math import *

WGS84_a  = 6378137.0
WGS84_b  = 6356752.314
WGS84_e_sq = 6.69438006676e-03
WGS84_ep_sq = 6.73949681994e-03

# Convert latitude, longitude, ellipsoidal height to ECEF Cartesian coords
# (lat and long in rad, height in m)
def ToCartes(GPSPt):
    CartesPt =  [0,0,0]
    LatRad = GPSPt[0]      # lat in rad
    LongRad = GPSPt[1]     # long in rad
    h = GPSPt[2]           # height in m         
	
    N = WGS84_a / sqrt(1 - WGS84_e_sq * sin(LatRad)*sin(LatRad))
    CartesPt[0] = (N+h)*cos(LatRad)*cos(LongRad)      # x in m
    CartesPt[1] = (N+h)*cos(LatRad)*sin(LongRad)      # y in m
    CartesPt[2] = (N*(1-WGS84_e_sq) + h)*sin(LatRad)  # z in m
    return CartesPt

def ToGeodetic(CartesPt):
    GeodPt = [0,0,0]
    x = CartesPt[0]
    y = CartesPt[1]
    z = CartesPt[2]
    
    p = sqrt(x*x + y*y)
    theta = atan(z*WGS84_a/(p*WGS84_b))
    phi = atan((z+WGS84_ep_sq*WGS84_b*sin(theta)*sin(theta)*sin(theta))/\
                     (p-WGS84_e_sq*WGS84_a*cos(theta)*cos(theta)*cos(theta)))
    N = WGS84_a / sqrt(1 - WGS84_e_sq * sin(phi)*sin(phi))
    GeodPt[0] = phi
    GeodPt[1] = atan2(y,x)
    print p/cos(phi)
    print N
    GeodPt[2] = p/cos(phi) - N
    return GeodPt    
    
# Calc Cartesian distance
def CartesDist(p1,p2):
    dist = (p2[0]-p1[0])*(p2[0]-p1[0]) + \
               (p2[1]-p1[1])*(p2[1]-p1[1]) + \
               (p2[2]-p1[2])*(p2[2]-p1[2])
    return sqrt(dist)

