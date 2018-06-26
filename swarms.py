"""A detailed work-through of the model for collisional evolution of irregular
satellite swarms presented in the paper Kennedy G.M., Wyatt M. C., 2011, Monthly
Notices of the Royal Astronomical Society, 412, 2137.

EVERY INPUT MUST BE IN SI UNITS

IN DEVELOPMENT"""

from numpy import pi, inf
import scipy.integrate as integrate

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
    Dmin: float; Dt: float; Dc: float; Ma: float
    Dmax: float; rho: float; qs: float; sigma0: float
    kg_val: float; ks_val: float; qg: float

    def __init__(self, Dmin, Dt, Dmax, Dc=None, Ma=None, sigma0=None, rho=1500, qs=1.9, qg=1.7):
        self.Dmin = Dmin; self.Dt = Dt; self.Dc = Dmax
        self.Dmax = Dmax; self.rho = rho; self.qs = qs
        self.qg = qg; self.kg_val = None; self.ks_val = None;
        if Dc is not None:
            self.Dc = Dc
        if Ma is not None:
            self.Ma = Ma
            self.init_from_mass()
        elif sigma0 is not None:
            self.sigma0 = sigma0
            self.init_from_area()
        else:
            raise ValueError("Mass and Area both unspecified")

    def init_from_mass(self):
        """initialize from the specified mass"""
        self.kg_val = (self.Ma*(6/(pi*self.rho))*(6 - 3*self.qg)
                        *self.Dc**(3*self.qg - 6))
        self.kg_to_ks()

    def init_from_area(self):
        """initialize from the specified area"""
        self.ks_val = ((4*(3*self.qs - 5)/pi)*self.sigma0
                        *self.Dmin**(3*self.qs - 5))
        self.ks_to_kg()

    def kg_to_ks(self):
        """description TBD"""
        self.ks_val = self.Dt**(3*self.qs - 3*self.qg)*self.kg_val

    def ks_to_kg(self):
        """description TBD"""
        self.kg_val = self.Dt**(3*self.qg - 3*self.qs)*self.ks_val

    def Mtot(self, Rcc0, t, tnleft, A, M_init):
        """The total mass of the swarm at some given time t."""
        if t <= tnleft:
            return (M_init)/(1 + Rcc0*t)
        else:
            return A*self.Dc**3

    def Atot(self):
        """Surface area of swarm. Output in AU"""
        return (pi/4)*265.3*self.ks_val*(self.Dmin**(5 - 3*self.qs)/(3*self.qs - 5))

    def Atot_int(self, dlow, dmid, dhigh):
        """Surface area of swarm. Using integration. Output in meters."""
        if (dlow is None) and (dmid is None):
            lower = (self.ks_val/(5 - 3*self.qs))*(self.Dt**(5 - 3*self.qs)
                                                - self.Dmin**(5 - 3*self.qs))

        elif dlow is None:
            if dmid < self.Dmin:
                raise ValueError("Invalid transition size input")
            lower = (self.ks_val/(5 - 3*self.qs))*(dmid**(5 - 3*self.qs)
                                                - self.Dmin**(5 - 3*self.qs))
        elif dmid is None:
            if dlow > self.Dt:
                raise ValueError("Invalid min size input")
            lower = (self.ks_val/(5 - 3*self.qs))*(self.Dt**(5 - 3*self.qs)
                                                - dlow**(5 - 3*self.qs))
        else:
            if dlow > dmid:
                raise ValueError("Invalid min and/or transition size input")
            lower = (self.ks_val/(5 - 3*self.qs))*(dmid**(5 - 3*self.qs)
                                                - dlow**(5 - 3*self.qs))

        if (dmid is None) and (dhigh is None):
            upper = (self.kg_val/(5 - 3*self.qg))*(self.Dc**(5 - 3*self.qg)
                                                    - self.Dt**(5 - 3*self.qg))

        elif dmid is None:
            if dhigh < self.Dt:
                raise ValueError("Invalid max size input")
            upper = (self.kg_val/(5 - 3*self.qg))*(digh**(5 - 3*self.qg)
                                                    - self.Dt**(5 - 3*self.qg))

        elif dhigh is None:
            if dmid > self.Dc:
                raise ValueError("Invalid transition size input")
            upper = (self.kg_val/(5 - 3*self.qg))*(self.Dc**(5 - 3*self.qg)
                                                    - dmid**(5 - 3*self.qg))
        else:
            if dmid > dhigh:
                raise ValueError("Invalid transition and/or max size input")
            upper = (self.kg_val/(5 - 3*self.qg))*(dhigh**(5 - 3*self.qg)
                                                    - dmid**(5 - 3*self.qg))

        #print("lower = {0:.5e}".format((pi/4)*lower))
        #print("upper = {0:.5e}".format((pi/4)*upper))
        return (pi/4)*(lower + upper)

    def Atot_mod(self):
        """A testing function for SI input"""
        part1 = (3/2)*(1/self.rho)*((6 - 3*self.qg)/(3*self.qs - 5))
        part2 = self.Dmin**(5 - 3*self.qs)*self.Dt**(3*self.qs - 3*self.qg)
        return part1*part2*self.Dc**(3*self.qg - 6)*self.Ma

    def Ntot_mod(self):
        a = self.kg_val*(2**(3*self.qg - 3) - 1)*self.Dc**(3 - 3*self.qg)
        return a/(3*self.qg - 3)

    def Ntot(self, dlow, dmid, dhigh):
        """Number of objects between input specification"""
        if (dlow is None) and (dmid is None):
            #print("ks = ", self.ks_val)
            lower = (self.ks_val/(3 - 3*self.qs))*(self.Dt**(3 - 3*self.qs)
                                                - self.Dmin**(3 - 3*self.qs))
        elif dlow is None:
            lower = (self.ks_val/(3 - 3*self.qs))*(dmid**(3 - 3*self.qs)
                                                - self.Dmin**(3 - 3*self.qs))
        elif dmid is None:
            lower = (self.ks_val/(3 - 3*self.qs))*(self.Dt**(3 - 3*self.qs)
                                                - dlow**(3 - 3*self.qs))
        else:
            lower = (self.ks_val/(3 - 3*self.qs))*(dmid**(3 - 3*self.qs)
                                                - dlow**(3 - 3*self.qs))

        if (dmid is None) and (dhigh is None):
            upper = (self.kg_val/(3 - 3*self.qg))*(self.Dc**(3 - 3*self.qg)
                                                - self.Dt**(3 - 3*self.qg))
        elif dmid is None:
            upper = (self.kg_val/(3 - 3*self.qg))*(dhigh**(3 - 3*self.qg)
                                                - self.Dt**(3 - 3*self.qg))
        elif dhigh is None:
            upper = (self.kg_val/(3 - 3*self.qg))*(self.Dc**(3 - 3*self.qg)
                                                - dmid**(3 - 3*self.qg))
        else:
            upper = (self.kg_val/(3 - 3*self.qg))*(dhigh**(3 - 3*self.qg)
                                                - dmid**(3 - 3*self.qg))

        print("kg_val = ", self.kg_val)
        print("lower = ", lower)
        print("upper = ", upper)
        return (lower + upper)

    def n(self, D):
        """TBA"""
        return self.ks_val*D**(2 - 3*self.qs) + self.kg_val*D**(2 - 3*self.qg)

class CollSwarm:
    """Represents the irregular satellite swarm of a given planet. You can
    compute the following information using the following methods:
        -computeAtot(): computes the size distribution's surface area

    === Attributes ===
    Ma: initial mass of the swarm [kg]
    Dt: transition particle size of the swarm [m]
    Dmax: maximum particle size of the swarm [m]
    L_s: luminosity of the primary [W] (3.828e26 W = 1 solar lum)
    M_s: mass of the primary [kg] (1.989e30 kg = 1 solar mass)
    M_pl: mass of the planet [kg] (5.972e24 kg = 1 earth mass)
    a_pl: semi-major axis of the planet [m] (1.496e11 m = 1 AU)
    R_pl: radius of the planet [m]
    eta: constant fraction of the Hill radius at which irregulars orbit
    rho: mass density of the swarm [kg/m^3]
    """
    Dt: float; Dmax: float; Nstr: float
    L_s: float; M_s: float; M_pl: float; Dc: float
    a_pl: float; R_pl: float; swarm: SizeDistribution
    Dmin: float; eta: float; fQ: float; f_vrel: float
    Rcc0: float; tnleft: float; M_init: float; correction: bool

    def __init__(self, M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr,
                rho=1500, fQ=5, f_vrel=4/pi, correction=True):
        self.Dt = Dt; self.Dmax = Dmax; self.Nstr = Nstr
        self.L_s = L_s; self.M_s = M_s; self.M_pl = M_pl; self.Dc = Dmax
        self.a_pl = a_pl; self.R_pl = R_pl; self.eta = eta; self.rho = rho
        self.Dmin = self.computeDmin()/1e6; self.fQ = fQ; self.f_vrel = f_vrel
        self.swarm = SizeDistribution(self.Dmin, self.Dt, self.Dmax, Ma=M0)
        self.Rcc0 = self.computeRCC(); self.tnleft = self.computetnleft()
        self.M_init = M0; self.correction = correction

    def computeDmin(self):
        """Compute the minimum sized object in the distribution."""
        a1 = (self.eta**0.5)*(self.L_s/3.828e26)
        a2 = self.rho*((self.M_pl/5.972e24)**(1/3))*((self.M_s/1.989e30)**(2/3))
        return 2e5*(a1/a2)

    def computeAtot(self, dlow=None, dmid=None, dhigh=None):
        """Compute the distribution's surface area."""
        return self.swarm.Atot_int(dlow, dmid, dhigh)

    def computeNtot(self, dlow=None, dmid=None, dhigh=None):
        """Return the distribution's number of particles."""
        return self.swarm.Ntot(dlow, dmid, dhigh)

    def computen(self, D):
        """TBA"""
        return self.swarm.n(D)

    def computeQd(self, D):
        """Compute the planetesimal strength of an object."""
        return 0.1*self.rho*(D/1000)**1.26/self.fQ

    def computeRCC(self):
        """Compute the rate of collision."""
        Qd = self.computeQd(self.Dmax)
        a = (self.swarm.Ma/5.972e24)*(self.M_s/1.989e30)**1.38*self.f_vrel**2.27
        b = (Qd**0.63*self.rho*(self.Dmax/1000)*(self.M_pl/5.972e24)**0.24*
            ((self.a_pl/1.496e11)*self.eta)**4.13)
        return 1.3e7*(a/b)

    def computeVrel(self):
        """Compute the mean relative velocity of collisions."""
        a = (4/pi)*516*(self.M_pl/5.972e24)**(1/3)*(self.M_s/1.989e30)**(1/6)
        return a/((self.eta*(self.a_pl/1.496e11))**0.5)

    def computeXc(self):
        """Computes the constant Xc for which objects of size XcDc can destroy
        objects of size Dc."""
        Qd = self.computeQd(self.Dmax)
        vrel = self.computeVrel()
        return (2*Qd/(vrel**2))**(1/3)

    def computetnleft(self):
        """Compute the time at which the first object is stranded."""
        Xc = self.computeXc()
        print("Xc = {0:.3e}".format(Xc))
        print("Dc b4 nval = {0:.3e}".format(self.Dc))
        nval = self.swarm.Ntot_mod()
        #nval = self.computeNtot(dlow=Xc*self.Dc, dmid=Xc*self.Dc, dhigh=self.Dc)
        print("nval = {0:.3e}".format(nval))
        return nval/(self.Rcc0*self.Nstr)

    def computeDc(self, t):
        """Compute the new Dc after stranding time."""
        if (t < self.tnleft) or (not self.correction):
            return self.Dmax
        else:
            a = (1 + 0.4*(t - self.tnleft)/self.tnleft)**(1.2)
            return self.Dmax/a

    def computeMtot(self, t):
        """Compute the total mass at a given time t."""
        y0 = self.swarm.Mtot(self.Rcc0, t, inf, 0, self.M_init)
        A = y0/self.Dmax**3
        Mt = self.swarm.Mtot(self.Rcc0, t, self.tnleft, A, self.M_init)

        return Mt

    def updateSwarm(self, t):
        """Description TBD"""
        Dct = self.computeDc(t)
        #print("Dct = {0:.3e}".format(Dct))
        Mt = self.computeMtot(t)
        #print("Mt = {0:.3e}".format(Mt/5.972e24))
        self.swarm = SizeDistribution(self.Dmin, self.Dt, self.Dmax, Dc=Dct, Ma=Mt)
        self.Dc = Dct
