from math import *
from numpy import *


def ADPfrac2cart(cell, fracADP):
    (a, b, c, alphadeg, betadeg, gammadeg) = cell

    alpha = alphadeg / 360 * 2 * pi
    beta = betadeg / 360 * 2 * pi
    gamma = gammadeg / 360 * 2 * pi
    V = a * b * c * sqrt(1 - (cos(alpha)) ** 2 - (cos(beta)) ** 2 - (cos(gamma)) ** 2 \
                         + 2 * cos(alpha) * cos(beta) * cos(gamma))
    astar = b * c * sin(alpha) / V
    bstar = a * c * sin(beta) / V
    cstar = a * b * sin(gamma) / V
    # v is the volume of the parallelepiped formed by the unit vectors a/|a|, b/|b|, b/|c|.
    v = V / (a * b * c)

    BETA = array([[a, 0, 0], [b * cos(gamma), b * sin(gamma), 0]
                     , [c * cos(beta), c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)
                      , c * v / sin(gamma)]])

    D = array([[astar, 0, 0], [0, bstar, 0], [0, 0, cstar]])
    cartADP = dot(transpose(BETA), dot(D, dot(fracADP, dot(D, BETA))))
    return cartADP
