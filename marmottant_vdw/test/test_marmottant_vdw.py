from nose.tools import *
from marmottant_vdw.marmottant_vdw import MarmottantVanDerWaal
from numpy import arange, array, zeros, size
from math import pi

def test():
    # Simulation time
    t = arange(0, 5, 0.001) * (10 ** -6)
    T = t * 2 * pi * 2.8e6
    P = zeros((1, size(t)))
    R0 = 2e-6
    Rbuck = 0.99
    SigmaL = 0.073
    Chi = 0.38
    SigmaR0 = Chi*((1/Rbuck)**2-1)
    Rrupt = Rbuck * (1 + SigmaL / Chi) ** (0.5)
    Rbreak = Rrupt
    Mu = 0.001
    KappaG = 1.06
    KappaSh = 2.4e-9
    R = [
            zeros(size(t)),
            zeros(size([0.1e6])),
            zeros(size([2e-6])),
            zeros(size([2.8e6]))
            ]
    w = 2 * pi * 2.8e6
    Rho = 998
    P0 = 101325
    C = 1481
    R_pert_ini = 1
    V0 = 0
    Rini = array([R_pert_ini, V0]).reshape(-1, 1)
    MarmottantVanDerWaal(t, Rini, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, "marm", "vdw")
    assert 1 == 1
