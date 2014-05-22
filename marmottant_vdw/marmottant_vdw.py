from numpy import matrix, zeros, array
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from math import pi

__version__ = '1.0'

class MarmottantVanDerWaal:
    def __init__(self, t, R, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, bubble_eq, gas):
        # Interpolate pulse
        self._t = t
        self._R = R
        self._T = T
        self._P = P
        self._w = w
        self._R0 = R0
        self._Rbuck = Rbuck
        self._Rrupt = Rrupt
        self._Rbreak = Rbreak
        self._KappaSh = KappaSh
        self._Chi = Chi
        self._SigmaR0 = SigmaR0
        self._Rho = Rho
        self._P0 = P0
        self._SigmaL = SigmaL
        self._C = C
        self._Mu = Mu
        self._KappaG = KappaG
        self._bubble_eq = bubble_eq
        self._gas = gas
        self._V0 = 0
        self._R_pert = 0

    def marmottant_vdw(self, R, tspan, t, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, bubble_eq, gas):
        Pac = interp1d(T, P)(t)

        Rn = matrix([zeros(len(Pac)), zeros(len(Pac))])

        if bubble_eq == "marm":
            if R[0] < Rbuck:
                SigmaR = 0
            elif R[0] >= Rbuck and R[0] < Rbreak:
                SigmaR = Chi * ((R[0] / Rbuck) ** 2 - 1)
            elif R[0] >= Rupt:
                SigmaR = SigmaL
            else:
                SigmaR = SigmaL
        elif bubble_eq == "dejong":
            SigmaR0 = SigmaL
            KappaSh = KappaSh * 16 * pi
            SigmaR = SigmaL + 2 * Chi * (R[0] * R0) * (1 / R0 - 1 / (R[0] * R0))
        elif bubble_eq == "free":
            KappaSh = 0
            SigmaR = SigmaL
            SigmaR0 = SigmaR

        if gas == "vdw":
            h = R0 / 5.6
            Rn[0] = R[1]

            Rn[1] = (1/(Rho * w ** 2 * R0 ** 2 * R[0])) * ((-3 / 2 * Rho * (R0 * w * R[1]) ** 2) + (P0 + 2 * SigmaR0 / R0) * (((R[0] * R0) ** 3 - h ** 3) / (R0 ** 3 - h ** 3)) ** (-KappaG) * (1 - 3 * KappaG / C * (R0 * w * R[1]) * (R[0] * R0) ** 3 / ((R[0] * R0) ** 3 - h ** 3)) - 2 * SigmaR / (R0 * R[0]) - 4 * Mu * w * R[1] / R[0] - 4 * KappaSh * w * R[1] / (R0 * R[0] ** 2) - (P0 + P0 * Pac))
        elif gas == "ideal":
            Rn[0] = R[1]

            Rn[1] = (1 / (Rho * w ** 2 * R0 ** 2 * R[0])) * ((-3 / 2 * Rho * (R0 * w * R[1]) ** 2) + (P0 + 2 * SigmaR0 / R0) * (1 / R[0]) ** (3 * KappaG) * (1 - 3 * KappaG / C * R0 * w * R[1]) - 2 * SigmaR / (R0 * R[0]) - 4 * Mu * w * R[1] / R[0] - 4 * KappaSh * w * R[1] / (R0 * R[0] ** 2) - (P0 + P0 * Pac))

        print("**" + str(Rn[0]))
        print("++" + str(Rn[1]))
        print(">>" + str(array(Rn.flatten(1))[0]))
        return array(Rn.flatten(1))[0]

    def solve(self):
        solution, yinfo = odeint(self.marmottant_vdw, self._R, self._t, args = (self._t, self._T, self._P, self._w, self._R0, self._Rbuck, self._Rrupt, self._Rbreak, self._KappaSh, self._Chi, self._SigmaR0, self._Rho, self._P0, self._SigmaL, self._C, self._Mu, self._KappaG, self._bubble_eq, self._gas), full_output = 1);
        return [solution, yinfo]
