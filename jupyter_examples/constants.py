SIG = 5.670367e-8 #Stefan-Boltzmann constant
G =  6.6743e-11 # Gravitational constant
C = 299792458 #speed of light
K_B = 1.38064852e-23 #Boltzmann's constant
H = 6.626070040e-34 #Plank's constant
MEARTH = 5.972e24 # kg
REARTH = 6.3781e6 # m
MJUP = 1.89813e27 # kg
RJUP = 71492000 # m
MMOON = 7.348e22 # kg
PC = 3.086e16 # m
AU = 1.496e11 # m
MSUN = 1.99e30 # kg
LSUN = 3.828e26 # watts
TSUN = 5800 # Kelvin
RSUN = 6.955e8  # m
MICRON = 1e-6 # m
KM = 1e3 # m
YEAR = 3.154e7 # seconds
JY = 1.e-26
G = 6.67e-11
GCC = 1000 # kg/m^3

# Some useful star types stats
STAR_TYPES = {
    "G2": {
        "L": 1.*LSUN,
        "M": 1.*MSUN,
        "T": 5780
    },
    "M0": {
        "L": 0.072*LSUN,
        "M": 0.60*MSUN,
        "T": 3800
    },
    "K5": {
        "L": 0.16*LSUN,
        "M": 0.69*MSUN,
        "T": 4410
    },
    "G5": {
        "L": 0.79*LSUN,
        "M": 0.93*MSUN,
        "T": 5610
    },
    "F5": {
        "L": 2.5*LSUN,
        "M": 1.3*MSUN,
        "T": 6540
    },
    "F0": {
        "L": 5.2*LSUN,
        "M": 1.4*MSUN,
        "T": 7420
    }
}
