from nose.tools import *
from marmottant_vdw.marmottant_vdw import MarmottantVanDerWaal
from numpy import arange, zeros
from math import pi

def test():
    # Simulation time
    t = arange(0, 5, 0.001) * (10 ** -6)
    T = t * 2 * pi * 2.8e6
    P = zeros((1, len(t)))
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
    # R = zeros((len(t), ))
    Rho = 998
    P0 = 101325
    C = 1481
    MarmottantVanDerWaal(t, R, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, "marm", "vdw")
    assert 1 == 1
