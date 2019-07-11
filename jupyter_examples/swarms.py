"""A detailed work-through of the model for collisional evolution of irregular
satellite swarms presented in the paper Kennedy G.M., Wyatt M. C., 2011, Monthly
Notices of the Royal Astronomical Society, 412, 2137.

EVERY INPUT MUST BE IN SI UNITS

IN DEVELOPMENT"""

from numpy import pi, inf, exp, zeros, log2, log10, sqrt
import scipy.integrate as integrate
from pread import BaraffeModelFixedTime, BaraffeModelFixedMass

SIG = 5.670367e-8 #Stefan-Boltzmann constant
C = 299792458 #speed of light
K_B = 1.38064852e-23 #Boltzmann's constant
H = 6.626070040e-34 #Plank's constant
MEARTH = 5.972e24 # kg
MSUN = 1.99e30 # kg
AU = 1.496e11 # m
LSUN = 3.828e26 # watts
MICRON = 1e-6 # m
KM = 1e3 # m
YEAR = 3.154e7 # seconds
JY = 1.e-26
G = 6.67e-11

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

    def __init__(self, Dmin, Dmax, Dc=None, M0=None, sigma0=None, Dt=100, rho=1500, qs=1.9, qg=1.7, Nstr=6):
        self.Dmin = Dmin; self.Dt = Dt; self.Dc = Dmax
        self.Dmax = Dmax; self.rho = rho; self.qs = qs
        self.qg = qg; self.kg_val = None; self.ks_val = None;
        self.Nkg = None; self.Nks = None;
        self.A1 = []; self.A2 = [];
        self.Nstr = Nstr;
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

    def compute_kg_from_Mtot(self, Mt):
        """description TBD"""
        return (Mt*(6/(pi*self.rho))*(6 - 3*self.qg)
                    * self.Dmax**(3*self.qg - 6))

    def compute_ks_from_kg(self):
        """description TBD"""
        return self.Dt**(3*self.qs - 3*self.qg)*self.Nkg

    def Mtot(self, dlow=None, dhigh=None):#, Nkg, Nks):
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

    def DMtot(self, Rcc0, t, tnleft, Xc_val, M_init, correction):
        """The total mass of the swarm at some given time t."""
        numerator = (3*self.qg - 3)*self.Nstr*pi*self.rho
        denominator = 6*(2**(3*self.qg - 3) - 1) * (6 - 3*self.qg)
        A = numerator/denominator
        if (t <= tnleft) or (not correction):
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
        return (pi/4)*(lower + upper)

    def Atot_mod(self):
        """A testing function for SI input"""
        part1 = (3/2)*(1/self.rho)*((6 - 3*self.qg)/(3*self.qs - 5))
        part2 = self.Dmin**(5 - 3*self.qs)*self.Dt**(3*self.qs - 3*self.qg)
        return part1*part2*self.Dc**(3*self.qg - 6)*self.M0

    def Ntot_mod(self):
        a = self.kg_val*(2**(3*self.qg - 3) - 1)*self.Dc**(3 - 3*self.qg)
        return a/(3*self.qg - 3)

    def Ntot(self, dlow, dhigh, verbose=None):
        """Number of objects between input specification"""
        if dlow is None:
            dlow = self.Dmin
            dmid= self.Dt
        elif dlow > self.Dt:
            dmid = dlow
        else:
            dmid = self.Dt
        if dhigh is None:
            dhigh = self.Dmax

        lower = (self.ks_val/(3 - 3*self.qs))*(self.Dt**(3 - 3*self.qs)
                                            - dlow**(3 - 3*self.qs))

        upper = 0
        str_upper = 0
        K_str = self.kg_val * self.Dc**(3 - 3*self.qg) #numerator / denominator
        str_upper = K_str * (log2(self.Dmax) - log2(self.Dc))
        upper = (self.kg_val/(3 - 3*self.qg))*(self.Dc**(3 - 3*self.qg)
                                            - dmid**(3 - 3*self.qg))

        if dlow > self.Dc:
            # compute K_str
            print("\t dlow > Dc")
            str_upper = K_str * (log2(self.Dmax) - log2(dlow))
            upper = 0
            lower = 0

        elif self.Dt < dlow <= self.Dc:
            print("\t dlow < Dc")
            upper = (self.kg_val/(3 - 3*self.qg))*(self.Dc**(3 - 3*self.qg)
                                                - dlow**(3 - 3*self.qg))
            str_upper = K_str * (log2(self.Dmax) - log2(self.Dc))
            lower = 0

        if (verbose is not None):
            from random import randint
            num = randint(0, 100)
            if num == 5:
                print("ks_val = ", self.ks_val)
                print("kg_val = ", self.kg_val)
                print("lower = ", lower)
                print("upper = ", upper)
                print("str upper = ", str_upper)
                print("qg = ", self.qg)
                print("k_str = ", K_str)
                print("Dmax = ", self.Dmax)
                print("Dc = ", self.Dc)
                print("dlow = ", dlow)
        if dlow > self.Dmax:
            return 0
        elif dlow > self.Dt:
            return (upper + str_upper)
        else:
            return (lower + upper + str_upper)

    def n(self, D):
        """TBA"""
        return self.ks_val*D**(2 - 3*self.qs) + self.kg_val*D**(2 - 3*self.qg)
    
def computeBnu(lamb, T):
    a = 2*H*(C/lamb)**3/C**2
    b = 1/(exp(H*(C/lamb)/(K_B*T)) - 1)
    return a*b
   
def computeFthermal(lamb, Across, T, d):
    Bnu = computeBnu(lamb, T)
    return Across*Bnu/d**2

def computeCRthermal(lamb, Tp, Tstar, Ap, Astar):
    Bnup = computeBnu(lamb, Tp)
    Bnustar = computeBnu(lamb, Tstar)
    return Ap/Astar*Bnup/Bnustar

def lum_to_temp(L, A):
    ''' Get temperature from luminosity and surface area

    === Arguments ===
    L: Bolometric luminosity (W)
    A: Surface area (m^2)
    '''
    return (L/A/SIG)**(1./4.)

def computeFscat(star, lamb, g, Q, A, dscat):
    ''' Calculate light scattered from scattering surface with area A to the observer at a distance star.d away from star

    === Attributes ===
    star: star whose light is being scattered
    lamb: wavelength of observation
    g: scattering phase function
    Q: Albedo at wavelength lamb
    A: Total scattering area
    dscat: Distance of scattering surface from the central star
    '''

    Fincident = star.F(lamb, dscat) # flux incident from star on scattering area
    return Fincident*g*Q*A/(pi*star.d**2)

def computeCRscat(g, Q, A, dscat):
    ''' Calculate contrast ratio in scattered light between object with area A and central star

    === Attributes ===
    g: scattering phase function
    Q: Albedo at wavelength lamb
    A: Total scattering area
    dscat: Distance of scattering surface from the central star
    '''

    return g*Q*A/(pi*dscat**2)

class Star():
    def __init__(self, L, M, T, d):
        self.L = L # Luminosity
        self.M = M # Mass
        self.T = T # Temperature
        self.d = d # distance from solar system
    
    def F(self, lamb, dist):
        # Flux from star at wavelength lamb at distance dist
        return computeFthermal(lamb, self.Across, self.T, dist)

    def mag(self, lamb, F0):
        F = self.F(lamb, self.d)
        return -2.5*log10(F/F0)

    def Imag(self): # Bessel (1979)
        lamb = 0.79*MICRON
        F0 = 2550*JY
        return self.mag(lamb, F0)
    
    def Hmag(self): # Campins et al. (1985)
        lamb = 1.6*MICRON
        F0 = 1080*JY
        return self.mag(lamb, F0)

    @property
    def A(self):        # 4piR**2
        return self.L/SIG/self.T**4
    
    @property
    def Across(self):
        return self.A/4.# piR**2
    
class Planet():
    def __init__(self, star, M, a, Q, R=None, Z='010', age=1.e10*YEAR):
        self.star = star
        self.M = M
        self.a = a
        self.Q = Q
        self.Z = Z
        self.age = age
        self.R = R if R is not None else self.RBaraffe(self.age) 
    
    @property
    def Lintrinsic(self):
        if self.M < 20*MEARTH: # minimum value in Baraffe. Bellow this Lintrinsic is negligible
            return 0.
        model = BaraffeModelFixedTime(self.Z, self.age)
        return model.L(self.M)
    
    @property
    def Lincident(self):
        return self.star.L/(4*self.a**2)*self.R**2*(1.-self.Q)

    @property
    def Ltot(self):
        return self.Lintrinsic + self.Lincident
    
    @property
    def T(self):
        return lum_to_temp(self.Ltot, self.A)

    @property
    def A(self):
        return 4.*pi*self.R**2
    
    @property
    def Across(self):
        return pi*self.R**2
   
    @property
    def RH(self):
        """ Hill radius """
        return self.a*(self.M/3./self.star.M)**(1./3.)

    def Fscat(self, lamb, g):
        # Calculate light scattered from planet to the observer
        return computeFscat(self.star, lamb, g, self.Q, pi*self.R**2, self.a)

    def CRscat(self, g):
        return computeCRscat(g, self.Q, pi*self.R**2, self.a)
        
    def RBaraffe(self, t):
        model = BaraffeModelFixedTime(self.Z, t)
        return model.R(self.M)
    
    def Fthermal(self, lamb):
        return computeFthermal(lamb, self.Across, self.T, self.star.d)
    
    def CRthermal(self, lamb):
        return computeCRthermal(lamb, self.T, self.star.T, self.A, self.star.A)
        
class CollSwarm:
    """Represents the irregular satellite swarm of a given planet. You can
    compute the following information using the following methods:
        -computeAtot(): computes the size distribution's surface area

    === Attributes ===
    M0: initial mass of the swarm [kg]
    Dt: transition particle size of the swarm [m]
    Dmax: maximum particle size of the swarm [m]
    eta: constant fraction of the Hill radius at which irregulars orbit
    rho: mass density of the swarm [kg/m^3]
    """
    Dt: float; Dmax: float; Nstr: float; Dc: float; swarm: SizeDistribution
    Dmin: float; eta: float; fQ: float; f_vrel: float; 
    Rcc0: float; tnleft: float; M_init: float; correction: bool; Dmin_min: float

    def __init__(self, star, planet, M0, Dt, Dmax, eta, Nstr, Q, rho=1500, fQ=5, f_vrel=4/pi, correction=True, alpha=1./1.2, Dmin_min = 1.65, age=0.):
        self.planet = planet; self.star = star
        self.Dt = Dt; self.Dmax = Dmax; self.Nstr = Nstr
        self.alpha = alpha; self.Dmin_min = Dmin_min
        self.age = age
        self.Dc = Dmax; self.eta = eta; self.rho = rho; self.Q = Q
        self.Dmin = self.computeDmin(self.Dmin_min); self.fQ = fQ; self.f_vrel = f_vrel
        self.swarm = SizeDistribution(self.Dmin, self.Dmax, M0=M0); self.M_init = M0
        self.Rcc0 = self.computeRCC(); self.tnleft = self.computetnleft()
        self.correction = correction
        self.updateSwarm(self.age)

    def computeDmin2(self, Dmin_min=None):
        """Compute the minimum sized object in the distribution."""
        if Dmin_min is None:
            Dmin_min = self.Dmin_min
        a1 = (self.eta**0.5)*(self.star.L/LSUN)
        a2 = self.rho*((self.planet.M/MEARTH)**(1/3))*((self.star.M/MSUN)**(2/3))
        return max(2e5*(a1/a2)*MICRON, Dmin_min)
    
    def computeDmin(self, Dmin_min=None):
        """Compute the minimum sized object in the distribution."""
        if Dmin_min is None:
            Dmin_min = self.Dmin_min
        a1 = (self.eta**0.5)*(self.star.L/LSUN)
        a2 = self.rho*((self.planet.M/MEARTH)**(1/3))*((self.star.M/MSUN)**(2/3))
        return max(2e5*(a1/a2), Dmin_min)*MICRON

    def computeAtot(self, dlow=None, dmid=None, dhigh=None, cap=False):
        """Compute the distribution's surface area."""
        if cap:
            R_H = self.planet.a * (self.planet.M / (3 * self.star.M))**(1./3.)
            return min(self.swarm.Atot(dlow, dhigh), pi * R_H**2)
        else:
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
        return 0.1*self.rho*(D/KM)**1.26/self.fQ

    @property
    def T(self):
        # Calculate equilibrium temperature of the planet from incident light
        return (self.star.L*(1.-self.Q)/(16*pi*SIG*self.planet.a**2))**(1./4.)
    
    def Fscat(self, lamb, g, planet=False, swarm=False, Fstar=None):
        return computeFscat(self.star, lamb, g, self.Q, self.computeAtot(), self.planet.a)

    def CRscat(self, g):
        return computeCRscat(g, self.Q, self.computeAtot(), self.planet.a)
        
    def Fthermal(self, lamb):
        return computeFthermal(lamb, self.computeAtot(), self.T, self.star.d)
    
    def CRthermal(self, lamb):
        return computeCRthermal(lamb, self.T, self.star.T, self.computeAtot(), self.star.A)

    def computeRCC(self):
        """Compute the rate of collision."""
        Qd = self.computeQd(self.Dc)
        a = (self.swarm.M0/MEARTH)*(self.star.M/MSUN)**1.38*self.f_vrel**2.27
        b = (Qd**0.63*self.rho*(self.Dmax/KM)*(self.planet.M/MEARTH)**0.24*
            ((self.planet.a/AU)*self.eta)**4.13)
        return 1.3e7*(a/b)/YEAR

    def computeVrel(self):
        """Compute the mean relative velocity of collisions."""
        return (4/pi)*sqrt(G*self.planet.M/(self.planet.RH*self.eta))
    
    def computeVrel2(self):
        """Compute the mean relative velocity of collisions."""
        a = (4/pi)*516*(self.planet.M/MEARTH)**(1/3)*(self.star.M/MSUN)**(1/6)
        return a/((self.eta*(self.planet.a/AU))**0.5)

    def computeXc(self):
        """Computes the constant Xc for which objects of size XcDc can destroy
        objects of size Dc."""
        Qd = self.computeQd(self.Dc)
        vrel = self.computeVrel()
        return (2*Qd/(vrel**2))**(1/3)

    def computeXcmax(self):
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
            a = (1 + 0.4*(t - self.tnleft)/self.tnleft)**(self.alpha)
            return self.Dmax/a

    def computeMtot(self, t):
        """Compute the total mass at a given time t."""
        Xc_max = self.computeXcmax()
        Mt = self.swarm.DMtot(self.Rcc0, t, self.tnleft, Xc_max, self.M_init, self.correction)
        return Mt

    def updateSwarm(self, t):
        """Description TBD"""
        Dct = self.computeDc(t)
        self.swarm.Dc = Dct
        self.Dc = Dct
        Mt = self.computeMtot(t)
        self.swarm = SizeDistribution(self.Dmin, self.Dmax, Dc=Dct, M0=Mt)

    def aopt(self, t):
        # t must be in years!
        t /= YEAR
        part1 = (self.star.M/MSUN)**0.33 * self.f_vrel**0.55
        part2 = (self.planet.M/MEARTH)**0.06 * self.computeQd(self.Dc)**0.15 * self.eta
        part3 = t * self.M_init/MEARTH / self.rho / (self.Dc/KM)
        return 50. * part1 / part2 * part3**0.24 * AU

