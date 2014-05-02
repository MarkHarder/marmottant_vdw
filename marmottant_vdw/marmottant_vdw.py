from numpy import array, size
from scipy.interpolate import interp1d

__version__ = '1.0'

class MarmottantVanDerWaal:
    def __init__(self, t, R, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, bubble_eq, gas):
        # Interpolate pulse
        self._Pac = interp1d(T, P)(t)
        self._Rn = array([[0.001 for i in range(size(self._Pac))], [0.001 for i in range(len(self._Pac))]])
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

    def solve(self):
        if self._bubble_eq == "marm":
            if self._R[0][0] < self._Rbuck:
                self._SigmaR = 0
            elif self._R[0][0] >= self._Rbuck and self._R[0] < self._Rbreak:
                self._SigmaR = self._Chi * ((self._R[0][0] / self._Rbuck) ** 2 - 1)
            elif self._R[0][0] >= self._Rupt:
                self._SigmaR = self._SigmaL
            else:
                self._SigmaR = self._SigmaL

        if self._gas == "vdw":
            h = self._R0 / 5.6
            self._Rn[0] = self._R[1]

            self._Rn[1] = (1/(self._Rho * self._w ** 2 * self._R0 ** 2 * self._R[0][0])) * ((-3 / 2 * self._Rho * (self._R0 * self._w * self._R[1][0]) ** 2) + (self._P0 + 2 * self._SigmaR0 / self._R0) * (((self._R[0][0] * self._R0) ** 3 - h ** 3) / (self._R0 ** 3 - h ** 3)) ** (-self._KappaG) * (1 - 3 * self._KappaG / self._C * (self._R0 * self._w * self._R[1][0]) * (self._R[0][0] * self._R0) ** 3 / ((self._R[0][0] * self._R0) ** 3 - h ** 3)) - 2 * self._SigmaR / (self._R0 * self._R[0][0]) - 4 * self._Mu * self._w * self._R[1][0] / self._R[0][0] - 4 * self._KappaSh * self._w * self._R[1][0] / (self._R0 * self._R[0][0] ** 2) - (self._P0 + self._P0 * self._Pac))

        self._Rn = self._Rn.flatten(1)
        print(self._Rn)
