from numpy import array, size
from scipy.interpolate import interp1d

__version__ = '1.0'

class MarmottantVanDerWaal:
    def __init__(self, t, R, T, P, w, R0, Rbuck, Rrupt, Rbreak, KappaSh, Chi, SigmaR0, Rho, P0, SigmaL, C, Mu, KappaG, bubble_eq, gas):
        # Interpolate pulse
        Pac = interp1d(T, P)(t)
        Rn = array([[0.001 for i in range(size(Pac))], [0.001 for i in range(len(Pac))]])

        if bubble_eq == "marm":
            if R[0][0] < Rbuck:
                SigmaR = 0
            elif R[0][0] >= Rbuck and R[0] < Rbreak:
                SigmaR = Chi * ((R[0][0] / Rbuck) ** 2 - 1)
            elif R[0][0] >= Rupt:
                SigmaR = SigmaL
            else:
                SigmaR = SigmaL

        if gas == "vdw":
            h = R0 / 5.6
            Rn[0] = R[1]

            Rn[1] = (1/(Rho * w ** 2 * R0 ** 2 * R[0][0])) * ((-3 / 2 * Rho * (R0 * w * R[1][0]) ** 2) + (P0 + 2 * SigmaR0 / R0) * (((R[0][0] * R0) ** 3 - h ** 3) / (R0 ** 3 - h ** 3)) ** (-KappaG) * (1 - 3 * KappaG / C * (R0 * w * R[1][0]) * (R[0][0] * R0) ** 3 / ((R[0][0] * R0) ** 3 - h ** 3)) - 2 * SigmaR / (R0 * R[0][0]) - 4 * Mu * w * R[1][0] / R[0][0] - 4 * KappaSh * w * R[1][0] / (R0 * R[0][0] ** 2) - (P0 + P0 * Pac))

        Rn = Rn.flatten(1)
        print(Rn)
