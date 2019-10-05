"""A detailed work-through of the model for collisional evolution of irregular
satellite swarms presented in the paper Kennedy G.M., Wyatt M. C., 2011, Monthly
Notices of the Royal Astronomical Society, 412, 2137.

EVERY INPUT MUST BE IN SI UNITS

IN DEVELOPMENT"""

from numpy import pi, inf, exp, zeros, log2, log10, sqrt
import scipy.integrate as integrate
from pread import BaraffeModelFixedTime, BaraffeModelFixedMass

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

    def __init__(self, Dmin, Dmax, Dc=None, M0=None, sigma0=None, Dt=100, rho=1000, qs=1.9, qg=1.7):
        """Dc defaults to Dmax. Dmax now just in there for forward compatbility. All masses, areas and numbers only computed up to size Dc"""
        self.Dmin = Dmin; self.Dt = Dt
        self.Dmax = Dmax; self.rho = rho; self.qs = qs
        self.qg = qg; self.kg_val = None; self.ks_val = None;
        self.Nkg = None; self.Nks = None;
        self.A1 = []; self.A2 = [];
        if Dc is None:
            Dc = Dmax
        self.Dc = Dc
        if self.Dmin > self.Dc:
            raise ValueError("Dmin must be smaller than Dc in swarms.SizeDistribution")
        if self.Dc > self.Dmax:
            raise ValueError("Dc must be <= Dmax in swarms.SizeDistribution")

        if M0 is not None:
            self.M0 = M0
            self.init_from_mass()
        elif sigma0 is not None:
            self.sigma0 = sigma0
            self.init_from_area()
        else:
            raise ValueError("Mass and Area both unspecified")

    def init_from_mass(self):
        """initialize from the specified mass."""
        if self.Dc < self.Dt:
            raise ValueError("swarms.sizeDistribution initialized with a maximum size Dc ({0:.1e}) < Dt ({0:.1e}). Implementation assumes Dc > Dt.".format(self.Dc, self.Dt))

        self.kg_val = (self.M0*(6/(pi*self.rho))*(6 - 3*self.qg)
                        *self.Dc**(3*self.qg - 6))
        self.kg_to_ks()

    def init_from_area(self):
        """initialize from the specified area"""
        if self.Dmin > self.Dt:
            raise ValueError("swarms.sizeDistribution initialized with a miminimum size Dmin > Dt. Implementation assumes Dmin < Dt.")
        self.ks_val = ((4*(3*self.qs - 5)/pi)*self.sigma0
                        *self.Dmin**(3*self.qs - 5))
        self.ks_to_kg()

    def kg_to_ks(self):
        """Calculate ks from kg from continuity"""
        self.ks_val = self.Dt**(3*self.qs - 3*self.qg)*self.kg_val

    def ks_to_kg(self):
        """Calculate kg from ks from continuity"""
        self.kg_val = self.Dt**(3*self.qg - 3*self.qs)*self.ks_val

    def Mtot(self, dlow=None, dhigh=None):
        """The total mass of the swarm. By default goes from Dmin to Dc"""
        if dlow is None:
            dlow = self.Dmin
        if dhigh is None:
            dhigh = self.Dc
        
        if dlow > dhigh:
            raise ValueError("dlow={0:.1e} must be < dhigh={1:.1e} in swarms.SizeDistribution.Mtot()".format(dlow, dhigh))
        
        if dhigh > self.Dc:
            raise ValueError("dhigh={0:.1e} must be <= Dc={1:.1e} in swarms.SizeDistribution.Mtot()".format(dhigh, Dc))

        lower = 0; upper = 0;

        if dlow < self.Dt:  # Use qs slope
            dmid = min(dhigh, self.Dt)
            lower = (self.ks_val/(6 - 3*self.qs))*(dmid**(6 - 3*self.qs)
                                            - dlow**(6 - 3*self.qs))

        if dhigh > self.Dt: # Use qg slope
            dmid = max(dlow, self.Dt)
            upper = (self.kg_val/(6 - 3*self.qg))*(dhigh**(6 - 3*self.qg)
                                            - dmid**(6 - 3*self.qg))

        return self.rho*pi*(lower + upper)/6

    def Atot(self, dlow=None, dhigh=None):
        """Surface area of swarm. Using integration. By default goes from Dmin to Dc."""
        if dlow is None:
            dlow = self.Dmin
        if dhigh is None:
            dhigh = self.Dc
        
        if dlow > dhigh:
            raise ValueError("dlow={0:.1e} must be < dhigh={1:.1e} in swarms.SizeDistribution.Atot()".format(dlow, dhigh))
        
        if dhigh > self.Dc:
            raise ValueError("dhigh={0:.1e} must be <= Dc={1:.1e} in swarms.SizeDistribution.Atot()".format(dhigh, Dc))

        lower = 0; upper = 0;

        if dlow < self.Dt:  # Use qs slope
            dmid = min(dhigh, self.Dt)
            lower = (self.ks_val/(5 - 3*self.qs))*(dmid**(5 - 3*self.qs)
                                            - dlow**(5 - 3*self.qs))

        if dhigh > self.Dt: # Use qg slope
            dmid = max(dlow, self.Dt)
            upper = (self.kg_val/(5 - 3*self.qg))*(dhigh**(5 - 3*self.qg)
                                            - dmid**(5 - 3*self.qg))

        return (pi/4)*(lower + upper)

    def Ntot(self, dlow=None, dhigh=None):
        """Number of objects between input specification. By default goes from Dmin to Dc."""
        if dlow is None:
            dlow = self.Dmin
        if dhigh is None:
            dhigh = self.Dc

        if dlow > dhigh:
            raise ValueError("dlow={0:.1e} must be < dhigh={1:.1e} in swarms.SizeDistribution.Ntot()".format(dlow, dhigh))

        if dhigh > self.Dc:
            raise ValueError("dhigh={0:.1e} must be <= Dc={1:.1e} in swarms.SizeDistribution.Ntot()".format(dhigh, Dc))

        lower = 0; upper = 0;
        
        if dlow < self.Dt:  # Use qs slope
            dmid = min(dhigh, self.Dt)
            lower = (self.ks_val/(3 - 3*self.qs))*(dmid**(3 - 3*self.qs)
                                            - dlow**(3 - 3*self.qs))

        if dhigh > self.Dt: # Use qg slope
            dmid = max(dlow, self.Dt)
            upper = (self.kg_val/(3 - 3*self.qg))*(self.Dc**(3 - 3*self.qg)
                                            - dmid**(3 - 3*self.qg))
        return lower + upper

    def n(self, D):
        """Differential size distribution, i.e. # of particles between size D and D+dD"""
        if D > self.Dc:
            raise ValueError("D={0:.1e} must be <= Dc={1:.1e} in swarms.SizeDistribution.n()".format(D, Dc))

        if D < self.Dt:
            return self.ks_val*D**(2 - 3*self.qs)
        else:
            return self.kg_val*D**(2 - 3*self.qg)
    
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

    Fincident = star.computeFthermal(lamb, dscat) # flux incident from star on scattering area
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
    
    def computeFthermal(self, lamb, dist):
        # Flux from star at wavelength lamb at distance dist
        return computeFthermal(lamb, self.Across, self.T, dist)

    def computemag(self, lamb, F0):
        F = self.computeFthermal(lamb, self.d)
        return -2.5*log10(F/F0)

    @property
    def Imag(self): # Bessel (1979)
        lamb = 0.79*MICRON
        F0 = 2550*JY
        return self.computemag(lamb, F0)
   
    @property
    def Hmag(self): # Campins et al. (1985)
        lamb = 1.6*MICRON
        F0 = 1080*JY
        return self.computemag(lamb, F0)

    @property
    def A(self):        # 4piR**2
        return self.L/SIG/self.T**4
    
    @property
    def Across(self):
        return self.A/4.# piR**2
    
class Planet():
    def __init__(self, star, M, a, Q=0.5, R=None, Z='002', age=1.e10*YEAR):
        self.star = star
        self.M = M
        self.a = a
        self.Q = Q
        self.Z = Z
        self.age = age
        self.R = R if R is not None else self.computeRBaraffe(self.age) 
    
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

    def computeFscat(self, lamb, g):
        # Calculate light scattered from planet to the observer
        return computeFscat(self.star, lamb, g, self.Q, pi*self.R**2, self.a)

    def computeCRscat(self, g):
        return computeCRscat(g, self.Q, pi*self.R**2, self.a)
        
    def computeRBaraffe(self, t):
        model = BaraffeModelFixedTime(self.Z, t)
        return model.R(self.M)
    
    def computeFthermal(self, lamb):
        return computeFthermal(lamb, self.Across, self.T, self.star.d)
    
    def computeCRthermal(self, lamb):
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
    stranding: whether to account for stranding and update Dc as a function of time
    """
    Dt: float; Dmax: float; Nstr: float; Dc: float; swarm: SizeDistribution
    Dmin: float; eta: float; fQ: float; f_vrel: float; 
    Rcc0: float; tnleft: float; M_init: float; stranding: bool; Dmin_min: float

    def __init__(self, star, planet, M0, Dmax=1e5, Dt=100, eta=0.3, Q=0.1, rho=1000, fQ=5, Nstr=6, f_vrel=4/pi, stranding=False, alpha=1./1.2, Dmin_min = 0., age=0., qs=1.9, qg=1.7, C1=0.1, C2=-1.9):
        self.planet = planet; self.star = star
        self.Dt = Dt; self.Dmax = Dmax; self.Nstr = Nstr
        self.alpha = alpha; self.Dmin_min = Dmin_min
        self.age = age; self.C1 = C1; self.C2 = C2
        self.eta = eta; self.rho = rho; self.Q = Q
        self.fQ = fQ; self.f_vrel = f_vrel
        self.M_init = M0
        self.stranding = stranding 
        self.qs = qs; self.qg = qg

        # these quantities are calculated but never change
        self.Dmin = self.computeDmin()
        self.swarm = SizeDistribution(Dmin=self.Dmin, Dmax=self.Dmax, M0=M0, Dt=Dt, rho=rho, qs=qs, qg=qg) # needed to init Rcc0
        self.Tcol0 = self.computeTcol(Mtot=M0, Dc=Dmax)
        self.tnleft = self.computetnleft(self.Dmax)
        self.updateSwarm(self.age)

    def computeDmin(self):
        """Compute the minimum sized object in the distribution."""
        # Assumes Qpr = 1
        num = 9.*self.star.L*self.eta**(1./2.)
        denom  = 8.*pi*G*self.star.M*C*self.rho # factor of 2 from Burns to get Diameter rather than radius
        return max(num/denom*(self.star.M/self.planet.M)**(1./3.), self.Dmin_min)

    def computeNtot(self, dlow=None, dhigh=None):
        """Return the distribution's number of particles."""
        return self.swarm.Ntot(dlow, dhigh)

    def computen(self, D):
        """TBA"""
        a = self.swarm.n(D)
        return a

    def computeQd(self, D):
        """Compute the strength of planetesimal of size D"""
        return 0.1*self.rho*(D/KM)**1.26/self.fQ
    
    def computeAtot(self, dlow=None, dhigh=None, cap=True):
        """Compute the distribution's surface area."""
        if dlow is None:
            dlow = self.computeDmin()
        if dhigh is None:
            dhigh = self.computeDc()
        A = self.swarm.Atot(dlow, dhigh)
        if cap:
            R_H = self.planet.a * (self.planet.M / (3 * self.star.M))**(1./3.)
            return min(A, pi * R_H**2)
        else:
            return A

    def computeTinc(self): # swarm does not have the contribution to temperature from internal energy planet has. Just calc from incident light
        return (self.star.L*(1.-self.Q)/(16*pi*SIG*self.planet.a**2))**(1./4.)
    
    def computeFscat(self, lamb, g, planet=False, swarm=False, Fstar=None):
        return computeFscat(self.star, lamb, g, self.Q, self.computeAtot(), self.planet.a)

    def computeCRscat(self, g):
        return computeCRscat(g, self.Q, self.computeAtot(), self.planet.a)
        
    def computeFthermal(self, lamb):
        return computeFthermal(lamb, self.computeAtot(), self.computeTinc(), self.star.d)
    
    def computeCRthermal(self, lamb):
        return computeCRthermal(lamb, self.computeTinc(), self.star.T, self.computeAtot(), self.star.A)

    def computeTcol(self, Mtot=None, Dc=None):
        return 1./self.computeRCC(Mtot, Dc)

    def computeVrel(self):
        """Compute the mean relative velocity of collisions."""
        return (4/pi)*sqrt(G*self.planet.M/(self.planet.RH*self.eta))
   
    def computeXc(self, Dc=None):
        """Computes the constant Xc for which objects of size XcDc can destroy
        objects of size Dc."""
        if Dc is None:
            Dc = self.computeDc()
        return (2*self.computeQd(Dc)/(self.computeVrel()**2))**(1/3)

    def computeDc(self):
        """Compute the new Dc if we want to correct for stranding and if after stranding time."""
        if (self.stranding == False) or (self.age <= self.tnleft): 
            return self.Dmax
        else:
            a = (1 + 0.4*(self.age - self.tnleft)/self.tnleft)**(self.alpha)
            return self.Dmax/a
    
    def computetnleft(self, Dc):
        Nhalf = self.swarm.Ntot(Dc/2., Dc)
        return Nhalf/(self.Nstr/self.Tcol0)

    def computeMtot(self, Dc=None):
        """Compute the total mass at a given time t."""
        if (self.stranding == False) or (self.age <= self.tnleft):
            return (self.M_init)/(1 + self.age/self.Tcol0)
        else:
            if Dc is None:
                Dc = self.computeDc()
            numerator = (3*self.qg - 3)*self.Nstr*pi*self.rho
            denominator = 6*(2**(3*self.qg - 3) - 1) * (6 - 3*self.qg)
            A = numerator/denominator
            return A*Dc**3

    def updateSwarm(self, t):
        """Description TBD"""
        self.age = t
        self.swarm = SizeDistribution(Dmin=self.computeDmin(), Dmax=self.Dmax, Dc=self.computeDc(), M0=self.computeMtot(), Dt=self.Dt, rho=self.rho)

    def computeaopt(self, t): # From our Eq 6 (factor to increase a by to reach Tcol = t)
        power = (1.- 2./3.*self.C2)/2.+3.# exact power \approx 4.13 in Kennedy paper. Use exact so that if we put planet at a=aopt, we actually get exact Tcol
        return self.planet.a*(t/self.Tcol0)**(1/power)
