"""A detailed work-through of the model for collisional evolution of irregular
satellite swarms presented in the paper Kennedy G.M., Wyatt M. C., 2011, Monthly
Notices of the Royal Astronomical Society, 412, 2137.

EVERY INPUT MUST BE IN SI UNITS

IN DEVELOPMENT"""

from numpy import pi, inf, exp, zeros
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
    Dmin: float; Dt: float; Dc: float; M0: float
    Dmax: float; rho: float; qs: float; sigma0: float
    kg_val: float; ks_val: float; qg: float

    def __init__(self, Dmin, Dmax, Dc=None, M0=None, sigma0=None, Dt=100, rho=1500, qs=1.9, qg=1.7):
        self.Dmin = Dmin; self.Dt = Dt; self.Dc = Dmax
        self.Dmax = Dmax; self.rho = rho; self.qs = qs
        self.qg = qg; self.kg_val = None; self.ks_val = None;
        if Dc is not None:
            self.Dc = Dc
        if M0 is not None:
            self.M0 = M0
            self.init_from_mass()
        elif sigma0 is not None:
            self.sigma0 = sigma0
            self.init_from_area()
        else:
            raise ValueError("Mass and Area both unspecified")

    def init_from_mass(self):
        """initialize from the specified mass"""
        self.kg_val = (self.M0*(6/(pi*self.rho))*(6 - 3*self.qg)
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

    def Mtot(self, dlow=None, dhigh=None):
        """The total mass of the swarm"""
        if dlow is None:
            dlow = self.Dmin
        if dhigh is None:
            dhigh = self.Dc

        lower = (self.ks_val/(6 - 3*self.qs))*(self.Dt**(6 - 3*self.qs)
                                            - dlow**(6 - 3*self.qs))

        upper = (self.kg_val/(6 - 3*self.qg))*(dhigh**(6 - 3*self.qg)
                                            - self.Dt**(6 - 3*self.qg))

        return self.rho*pi*(lower + upper)/6

    def DMtot(self, Rcc0, t, tnleft, A, M_init):
        """The total mass of the swarm at some given time t."""
        if t <= tnleft:
            return (M_init)/(1 + Rcc0*t)
        else:
            return A*self.Dc**3

    def Atot_mod2(self):
        """Surface area of swarm. Output in AU. Testing function"""
        return (pi/4)*265.3*self.ks_val*(self.Dmin**(5 - 3*self.qs)/(3*self.qs - 5))

    def Atot(self, dlow=None, dhigh=None):
        """Surface area of swarm. Using integration. Output in meters."""
        if dlow is None:
            dlow = self.Dmin
        if dhigh is None:
            dhigh = self.Dc

        lower = (self.ks_val/(5 - 3*self.qs))*(self.Dt**(5 - 3*self.qs)
                                            - dlow**(5 - 3*self.qs))

        upper = (self.kg_val/(5 - 3*self.qg))*(dhigh**(5 - 3*self.qg)
                                            - self.Dt**(5 - 3*self.qg))

        #print("lower = {0:.5e}".format((pi/4)*lower))
        #print("upper = {0:.5e}".format((pi/4)*upper))
        return (pi/4)*(lower + upper)

    def Atot_mod(self):
        """A testing function for SI input"""
        part1 = (3/2)*(1/self.rho)*((6 - 3*self.qg)/(3*self.qs - 5))
        part2 = self.Dmin**(5 - 3*self.qs)*self.Dt**(3*self.qs - 3*self.qg)
        return part1*part2*self.Dc**(3*self.qg - 6)*self.M0

    def Ntot_mod(self):
        a = self.kg_val*(2**(3*self.qg - 3) - 1)*self.Dc**(3 - 3*self.qg)
        return a/(3*self.qg - 3)

    def Ntot(self, dlow, dhigh):
        """Number of objects between input specification"""
        if dlow is None:
            dlow = self.Dmin
        if dhigh is None:
            dhigh = self.Dc

        lower = (self.ks_val/(3 - 3*self.qs))*(self.Dt**(3 - 3*self.qs)
                                            - dlow**(3 - 3*self.qs))

        upper = (self.kg_val/(3 - 3*self.qg))*(dhigh**(3 - 3*self.qg)
                                            - self.Dt**(3 - 3*self.qg))

        from random import randint
        num = randint(0, 100)
        if num == 5:
            print("ks_val = ", self.ks_val)
            print("kg_val = ", self.kg_val)
            print("lower = ", lower)
            print("upper = ", upper)
            print("qg = ", self.qg)
        if dlow > self.Dt:
            return upper
        else:
            return (lower + upper)

    def n(self, D):
        """TBA"""
        return self.ks_val*D**(2 - 3*self.qs) + self.kg_val*D**(2 - 3*self.qg)

class CollSwarm:
    """Represents the irregular satellite swarm of a given planet. You can
    compute the following information using the following methods:
        -computeAtot(): computes the size distribution's surface area

    === Attributes ===
    M0: initial mass of the swarm [kg]
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
    Dmin: float; eta: float; fQ: float; f_vrel: float; d_pl: float
    Rcc0: float; tnleft: float; M_init: float; correction: bool

    def __init__(self, M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr, d_pl,
                rho=1500, fQ=5, f_vrel=4/pi, correction=True):
        self.Dt = Dt; self.Dmax = Dmax; self.Nstr = Nstr
        self.L_s = L_s; self.M_s = M_s; self.M_pl = M_pl; self.Dc = Dmax
        self.a_pl = a_pl; self.R_pl = R_pl; self.eta = eta; self.rho = rho
        self.Dmin = self.computeDmin()/1e6; self.fQ = fQ; self.f_vrel = f_vrel
        self.swarm = SizeDistribution(self.Dmin, self.Dmax, M0=M0)
        self.M_init = M0; self.d_pl = d_pl
        self.Rcc0 = self.computeRCC(); self.tnleft = self.computetnleft()
        self.correction = correction

    def computeDmin(self):
        """Compute the minimum sized object in the distribution."""
        a1 = (self.eta**0.5)*(self.L_s/3.828e26)
        a2 = self.rho*((self.M_pl/5.972e24)**(1/3))*((self.M_s/1.989e30)**(2/3))
        return 2e5*(a1/a2)

    def computeAtot(self, dlow=None, dmid=None, dhigh=None):
        """Compute the distribution's surface area."""
        return self.swarm.Atot(dlow, dhigh)

    def computeNtot(self, dlow=None, dhigh=None):
        """Return the distribution's number of particles."""
        return self.swarm.Ntot(dlow, dhigh)

    def computen(self, D):
        """TBA"""
        a = self.swarm.n(D)
        return a

    def computeQd(self, D):
        """Compute the planetesimal strength of an object."""
        return 0.1*self.rho*(D/1000)**1.26/self.fQ

    def computeT(self, L, dist):
        sig = 5.670367e-8 #Stefan-Boltzmann constant
        return (L/(16*pi*sig*(dist)**2))**(1./4.)

    def computeBmu(self, lamb, T):
        c = 299792458 #speed of light
        k_B = 1.38064852e-23 #Boltzmann's constant
        h = 6.626070040e-34 #Plank's constant
        a = 2*h*(c/lamb)**3/c**2
        b = 1/(exp(h*(c/lamb)/(k_B*T)) - 1)
        return a*b

    def computeFth(self, lamb, planet=False, swarm=False):
        if planet:
            T = self.computeT(self.L_s, self.a_pl)
            Bmu = self.computeBmu(lamb, T)
            Fth = Bmu*pi*(self.R_pl/(self.d_pl))**2
            return Fth
        if swarm:
            T = 278.3*(self.L_s/3.828e26)**(1./4.)/(6.68459e-12*self.a_pl)**0.5
            Bmu = self.computeBmu(lamb, T)
            A = self.computeAtot()
            Fth = zeros(len(lamb))
            for i in range(len(lamb)):
                if lamb[i] >= 0.00021:
                    Fth[i] = (0.00021/lamb[i])*(Bmu[i]*A)/(self.d_pl**2)
                else:
                    Fth[i] = Bmu[i]*A/self.d_pl**2
            return Fth

    def computeFstar(self, Bmu, T):
        sig = 5.670367e-8 #Stefan-Boltzmann constant
        a = self.L_s*Bmu
        b = 4*sig*(T**4)*((self.a_pl)**2)
        return a/b

    def stellarTemp(self):
        Ms = self.M_s/1.989e30
        if Ms >= 16:
            return 3e4
        elif 2.1 <= Ms < 16:
            return 2e4
        elif 1.4 <= Ms < 2.1:
            return 9e3
        elif 1.04 <= Ms < 1.4:
            return 7e3
        elif 0.8 <= Ms < 1.04:
            return 5.6e3
        elif 0.45 <= Ms < 0.8:
            return 4.5e3
        elif 0.08 <= Ms < 0.45:
            return 3e3

    def computeFs(self, lamb, g, Q, planet=False, swarm=False):
        T = self.stellarTemp()
        Bmu = self.computeBmu(lamb, T)
        Fstar = self.computeFstar(Bmu, T)

        if planet:
            a = Fstar*self.R_pl**2*g*Q
            b = (self.d_pl)**2
            return a/b
        if swarm:
            a = Fstar*self.computeAtot()*g*Q
            b = pi*(self.d_pl)**2
            return a/b

    def computeRCC(self):
        """Compute the rate of collision."""
        Qd = self.computeQd(self.Dmax)
        a = (self.swarm.M0/5.972e24)*(self.M_s/1.989e30)**1.38*self.f_vrel**2.27
        b = (Qd**0.63*self.rho*(self.Dmax/1000)*(self.M_pl/5.972e24)**0.24*
            ((self.a_pl/1.496e11)*self.eta)**4.13)
        return 1.3e7*(a/b)

    def computeRCC5(self):
        """Compute the rate of collision."""
        Qd = self.computeQd(self.Dmax)
        a = (self.swarm.M0)*(self.M_s)**1.38*self.f_vrel**2.27
        b = (Qd**0.63*self.rho*(self.Dmax)*(self.M_pl)**0.24*
            ((self.a_pl)*self.eta)**4.13)
        return 39.16*(a/b)

    def computeRCC2(self):
        """Compute RCC using equation 5 from Kennedy instead of equation 7."""
        vrel = self.computeVrel()
        Xc = self.computeXc()
        V = 4*pi*(self.eta**2)*(self.eta/2)*((self.R_pl/1.496e11)**3)*0.866
        a = (6 - 3*self.swarm.qg)/(3*self.swarm.qs - 5)
        b = (vrel*0.1*(Xc**-1.9)*self.M_init)/(self.rho*(self.Dc/1000)*V)
        return 8.4e-5*a*b

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
        nval = self.swarm.Ntot_mod()
        return nval/(self.Rcc0*self.Nstr)

    def computeDc(self, t):
        """Compute the new Dc after stranding time."""
        if (t < self.tnleft) or (not self.correction):
            return self.Dmax
        else:
            a = (1 + 0.4*(t - self.tnleft)/self.tnleft)**(1/1.2)
            return self.Dmax/a

    def computeMtot(self, t):
        """Compute the total mass at a given time t."""
        y0 = self.swarm.DMtot(self.Rcc0, t, inf, 0, self.M_init)
        A = y0/self.Dmax**3
        Mt = self.swarm.DMtot(self.Rcc0, t, self.tnleft, A, self.M_init)

        return Mt

    def updateSwarm(self, t):
        """Description TBD"""
        Dct = self.computeDc(t)
        #print("Dct = {0:.3e}".format(Dct))
        Mt = self.computeMtot(t)
        #print("Mt = {0:.3e}".format(Mt/5.972e24))
        self.swarm = SizeDistribution(self.Dmin, self.Dmax, Dc=Dct, M0=Mt)
        self.Dc = Dct
