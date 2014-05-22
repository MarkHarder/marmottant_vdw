from nose.tools import *
from marmottant_vdw.marmottant_vdw import MarmottantVanDerWaal
from numpy import arange, array, zeros, size, sin
from math import pi
import matplotlib as mpl
import matplotlib.pyplot as plt

def test():
    # Simulation time
    t = arange(0, 5, 0.001) * (10 ** -6)
    T = t * 2 * pi * 2.8e6
    P = zeros((1, size(t)))
    R0 = 2*10e-6
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
    R_pert_ini = 1.001
    V0 = 0
    Rini = array([R_pert_ini, V0])

    pulse_type = "sin"
    Cycles = 3
    Period = 2 * pi
    T0 = Period * Cycles
    PRP_vect = 0.1e6

    if pulse_type == "sin":
        pulse_window = array(zeros(len(T)))
        for i, j in enumerate(T):
            if j > 0 and j < T0:
                pulse_window[i] = 1
            else:
                pulse_window[i] = 0

    if pulse_type == "kzk":
        pass
    else:
        PRP = PRP_vect
        A = PRP / P0
        P = -A * sin(w * t) * pulse_window

    solver = MarmottantVanDerWaal(t, Rini, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, "free", "ideal")
    solution, yinfo = solver.solve()
    print(solution[:,0])
    mpl.rc('font', family='serif', size=16)
    mpl.rc('xtick',labelsize='small')
    mpl.rc('ytick',labelsize='small')
    mpl.rc('legend',fontsize='small')
    plt.figure(1)
    plt.plot(t,solution[:,0],'bx')
    plt.legend(('Numerical','Exact'))
    plt.xlabel(r'Time [s]')
    plt.ylabel(r'Position [m]')
    plt.show()
    assert 1 == 1
