"""A detailed work-through of the model for collisional evolution of irregular
satellite swarms presented in the paper Kennedy G.M., Wyatt M. C., 2011, Monthly
Notices of the Royal Astronomical Society, 412, 2137.

EVERY INPUT MUST BE IN SI UNITS

IN DEVELOPMENT"""

from numpy import pi

class SizeDistribution:
    """An object that encapsulate the size distribution of the collision swarm.
    It includes information such as the total mass of the swarm, the surface
    area distribution, and the number of particles distribution at the size
    interval [X_c*D_c, D_c] and [Dmin, Dt].

    The size distribution object is not intended to be initialized on its own,
    but used with the CollSwarm object.

    === Attributes ===
    Dmin: minimum particle size of the swarm [m]
    Dt: transition particle size of the swarm [m]
    Dc: size of the largest non-stranded object [m]
    Dmax: maximum particle size of the swarm [m]
    rho: mass density of the swarm [kg/m^3]
    qs: size distribution slope for objects of size < Dt
    qg: size distribution slope for objects of size > Dt
    """
    Dmin: float; Dt: float; Dc: float
    Dmax: float; rho: float; qs: float
    kg_val: float; ks_val: float; qg: float

    def __init__(self, Dmin, Dt, Dmax, Dc=0, rho=1500, qs=1.9, qg=1.7):
        self.Dmin = Dmin; self.Dt = Dt; self.Dc = Dmax
        self.Dmax = Dmax; self.rho = rho; self.qs = qs
        self.qg = qg; self.kg_val = 0; self.ks_val = 0;

    def kg(self, M0):
        """Size distribution proportionality constant for objects of size > Dt"""
        self.kg_val = ((M0/5.97219e24)*(6/(pi*self.rho))*(6 - 3*self.qg)
                        *self.Dc**(3*self.qg - 6))

    def ks(self):
        """Size distribution proportionality constant for objects of size < Dt"""
        self.ks_val = self.Dt**(3*self.qs - 3*self.qg)*self.kg_val

    def Mtot(self, M0, Rcc0, t, tnleft, A, Dc):
        """The total mass of the swarm at some given time t."""
        if t <= tnleft:
            return (M0/5.97219e24)/(1 + Rcc0*t)
        else:
            return A*Dccc**3

    def Atot(self):
        """Surface area of swarm. Output in AU"""
        return (pi/4)*265.3*self.ks_val*(self.Dmin**(5 - 3*self.qs)/(3*self.qs - 5))

    def Ntot(self):
        """Number of objects between Dmin and Dt"""
        return -self.ks_val*self.Dmin**(3 - 3*self.qs)/(3 - 3*self.qs)

    def n(self, Xc):
        """Number of objects between XcDc and Dc"""
        return self.__kg*self.Dc**(3 - 3*self.qg)/(
                3 - 3*self.qg)*(1 - Xc**(3 - 3*self.qg))

class CollSwarm:
    """Represents the irregular satellite swarm of a given planet. You can
    compute the following information using the following methods:
        -computeAtot(): computes the size distribution's surface area


    === Attributes ===
    Ma: initial mass of the swarm [kg]
    Dt: transition particle size of the swarm [m]
    Dmax: maximum particle size of the swarm [m]
    L_s: luminosity of the primary [solar lum] --TO CHANGE
    M_s: mass of the primary [solar mass] --TO CHANGE
    M_pl: mass of the planet [earth masses] --TO CHANGE
    a_pl: semi-major axis of the planet [AU] --TO CHANGE
    R_pl: radius of the planet [m]
    eta: constant fraction of the Hill radius at which irregulars orbit
    rho: mass density of the swarm [kg/m^3]
    """
    Ma: float; Dt: float; Dmax: float
    L_s: float; M_s: float; M_pl: float
    a_pl: float; R_pl: float; swarm: SizeDistribution
    Dmin: float; eta: float;

    def __init__(self, M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, rho=1500):
        self.Ma = M0; self.Dt = Dt; self.Dmax = Dmax;
        self.L_s = L_s; self.M_s = M_s; self.M_pl = M_pl;
        self.a_pl = a_pl; self.R_pl = R_pl; self.eta = eta; self.rho = rho
        self.Dmin = self.computeDmin()/1e6;
        self.swarm = SizeDistribution(self.Dmin, self.Dt, self.Dmax)

    def computeDmin(self):
        """ATM: inputs are in solar masses, luminosity and earth-masses."""
        a1 = (self.eta**0.5)*self.L_s
        a2 = self.rho*(self.M_pl**(1/3))*(self.M_s**(2/3))
        return 2e5*(a1/a2)

    def computeAtot(self):
        """description TBD"""
        self.swarm.kg(self.Ma)
        self.swarm.ks()
        return self.swarm.Atot()

#main()
def main():
    jupiter = CollSwarm(7.37307e19, 100, 150000, 1, 1, 317.46, 5.2, 6.9911e7, 0.4)
    area = jupiter.computeAtot()
    print('area = {0:.3e}'.format(area))

if __name__ == '__main__':
    main()
