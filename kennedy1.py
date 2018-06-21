from numpy import pi, linspace, exp, logspace, concatenate, transpose, array
import matplotlib.pyplot as plt
from scipy import integrate

qg = 1.7
qs = 1.9
c = 299792458 #speed of light
k_B = 1.38064852e-23 #Boltzmann's constant
h = 6.626070040e-34 #Plank's constant
sig = 5.670367e-8 #Stefan-Boltzmann constant
G = 6.67408e-11 #gravitational constant

def Mtot(rho, kg_val, Dt, Dc):
    a = (Dc**(6 - 3*qg))/(6 - 3*qg)
    return (pi/6)*rho*kg_val*a

def n1(D, K):
    return 5.97e24*K*D**(3 - 3*qg)/(3 - 3*qg)*(1 - 0.5**(3 - 3*qg))

def n_ks(K, D):
    return 5.97e24*K*D**(2 - 3*qs)

def N_SI(K, D):
    return -K*D**(3 - 3*qs)/(3 - 3*qs)

def n_kg(K, D):
    return 5.97e24*K*D**(2 - 3*qg)

#def ks(kg_val, Dt):
#    return kg_val*Dt**(3*qs - 3*qg)

def ks1(Dt, kg_val):
    return 265.3*Dt**(3*qs - 3*qg)*kg_val

def ks1SI(Dt, kg_val):
    return Dt**(3*qs - 3*qg)*kg_val

def Atot(ks_val, Dmin):
    return (pi/4)*ks_val*(Dmin**(5 - 3*qs)/(3*qs - 5))

def kg(Mtot0, Dc, rho):
    return Mtot0*(6/(pi*rho))*(6 - 3*qg)*Dc**(3*qg - 6)

def D_min(eta, L_s, rho, M_pl, M_s):
    """Minimum object size in swarm.
    L_s in solar lum
    rho in kg/m^3
    M_pl in Earth mass
    M_s in Solar mass"""
    #D_min in micro meter
    a1 = (eta**0.5)*L_s
    a2 = rho*(M_pl**(1/3))*(M_s**(2/3))
    return 2e5*(a1/a2)

def num_strand(K, D):
    a = K*(2**(3*qg - 3) - 1)*(1e3*D)**(3 - 3*qg)
    return a/(3*qg - 3)

def v_rel(M_pl, M_s, eta, a_pl):
    a = (4/pi)*516*M_pl**(1/3)*M_s**(1/6)
    return a/((eta*a_pl)**0.5)

def Vol_satellite(eta, deta, rH):
    return 4*pi*eta**2*deta*rH**3

def Rh(M_p, M_s, a):
    return a*(M_p/(3*M_s))**(1/3)

def Q_d(D, rho):
    a = 0.1*rho*(D/1000)**1.26
    return a/5

def Xc(Qd, vrel):
    return (2*Qd/(vrel**2))**(1/3)

def Rcc(vrel, C1, X_c, C2, Mtot0, rho, Dc, V):
    a = (6 - 3*qg)/(3*qg - 5)
    b = (vrel*C1*X_c**C2*Mtot0)/(rho*(Dc/1000)*V)
    return 8.4e-5*a*b

def Rcc1(Mtot0, M_s, f_vrel, Qd, rho, Dc, M_pl, a_pl, eta):
    a = Mtot0*M_s**1.38*f_vrel**2.27
    b = Qd**0.63*rho*(Dc/1000)*M_pl**0.24*(a_pl*eta)**4.13
    return 1.3e7*(a/b)

def Dct(Dmax, tnleft, t):
    if t < tnleft:
        return Dmax
    else:
        a = (1 + 0.4*(t - tnleft)/tnleft)**(1/1.2)
        return Dmax/a

def Mtott(Mtot0, Rcc0, t, tnleft, A, Dccc):
    #print('Mtot0 = {0:.3e}'.format(Mtot0))
    #print('Rcc0 = {0:.3e}'.format(Rcc0))
    #print('t = {0:.3e}'.format(t))
    if t <= tnleft:
        a = Mtot0/(1 + Rcc0*t)
    else:
        a = A*Dccc**3
    #print('result = {0:.3e}'.format(a))
    return a

def t_col(eta, a_pl, M_pl, M_s, atot):
    a = (eta*a_pl)**(7/2)*M_pl**(2/3)
    b = M_s**(7/6)*atot
    return 1e-5*(a/b)

def Xpr1(rho, Dmin, area, M_s, a_pl, eta, M_pl, L_s):
    a = rho*Dmin*area*M_s**(7/6)
    b = a_pl**(3/2)*eta**(7/2)*M_pl**(2/3)*L_s
    return 4e4*(a/b)

def Xpr2(area, M_s, a_pl, eta, M_pl):
    a = area*(M_s**0.5)
    b = a_pl**(3/2)*(eta**3)*M_pl
    return 8e9*(a/b)

#Spectral intensity
def B_mu(lamb, T):
    #const1 = 2/(lamb**2)
    #return const1*k_B*T
    a = 2*h*(c/lamb)**3/c**2
    b = 1/(exp(h*(c/lamb)/(k_B*T)) - 1)
    return a*b

#Thermal radiation
def F_th(spec, A, d):
    return (spec*A)/(d**2)

#Effective temperature
def T(L, a_pl):
    #return (L/(A*sig))**0.25
    return (L/(16*pi*sig*(1.5e11*a_pl)**2))**(1/4)

def F_th_planet(lamb, L_s, R_pl, d_pl):
    T_planet = T(3.827e26*L_s, d_pl)
    Bmu_planet = B_mu(lamb, T_planet)
    #Fth_planet = F_th(Bmu_planet, 4*pi*R_pl**2, d_pl)
    #Fth_planet = pi*R_pl**2*(Bmu_planet/(4*pi*(1.5e11*d_pl)**2))*(1/(4*pi*(1.5e11*d_pl)**2))
    Fth_planet = Bmu_planet*pi*(R_pl/(1.5e11*d_pl))**2

    return T_planet, Bmu_planet, Fth_planet

def F_th_dust(lamb, L_s, A, d_pl):
    T_dust = T(3.827e26*L_s, d_pl)
    Bmu_dust = B_mu(lamb, T_dust)
    Fth_dust = (0.00021/lamb)*F_th(Bmu_dust, A, d_pl)

    return T_dust, Bmu_dust, Fth_dust

def F_s(L_s, Bnu, T_s, a_pl):
    a = 3.827e26*L_s*Bnu
    b = 4*sig*(T_s**4)*((1.5e11*a_pl)**2)
    return a/b

def F_scat(Fs, R_pl, g, Q, d_pl):
    a = Fs*R_pl**2*g*Q
    b = (1.5e11*d_pl)**2
    return a/b

def F_scat_dust(Fs, A, g, Q, d_pl):
    a = Fs*A*g*Q
    b = pi*(d_pl)**2
    return a/b

def stellarTemp(M_s):
    if M_s >= 16:
        return 3e4
    elif 2.1 <= M_s < 16:
        return 2e4
    elif 1.4 <= M_s < 2.1:
        return 9e3
    elif 1.04 <= M_s < 1.4:
        return 7e3
    elif 0.8 <= M_s < 1.04:
        return 5.6e3
    elif 0.45 <= M_s < 0.8:
        return 4.5e3
    elif 0.08 <= M_s < 0.45:
        return 3e3

def g(y, t, params):
    const, tnleftt, Dmaxx = params[0], params[1], params[2]
    deriv = -(const/(Dmaxx**1.7938))*(1 + 0.4*(t - tnleftt)/tnleftt)**(1.7938*(1/1.2))*(y**2)
    return deriv

def main(Mtot0, Dt, Dc, rho, eta, L_s, M_pl, M_s, a_pl, R_pl):
    Dmax = Dc
    Dmin = D_min(eta, L_s, rho, M_pl, M_s)
    kg_val = kg(Mtot0, Dc, rho)
    ks_val = ks1(Dt, kg_val)
    #ks_val = ks(kg_val, Dt)
    area = Atot(ks_val, (Dmin/1e6))

    print('mtot0 = ', Mtot0)
    print('rho = {0:.3e}'.format(rho))

    Qd = Q_d(Dc, rho)
    vrel = v_rel(M_pl, M_s, eta, a_pl)
    X_c = Xc(Qd, vrel)
    rH = Rh(M_pl, 334672.021419*M_s, a_pl)
    V = Vol_satellite(eta, eta/2, rH)
    #R_cc = Rcc(vrel, 0.1, X_c, -1.9, Mtot0, rho, Dc, V)
    n = n1(Dc, kg_val)
    R_cc1 = Rcc1(Mtot0, M_s, 4/pi, Qd, rho, Dc, M_pl, a_pl, eta)
    tnleft = (n/(R_cc1*6))/1e6

    M_tott = Mtott(Mtot0, R_cc1, 4.5e9, 1e111, 1e10, 0)
    Dc_t = Dct(Dc, 1e6*tnleft, 4.5e9)
    kg_valt = kg(M_tott, Dc_t, rho)
    ks_valt = ks1(kg_valt, Dt)
    areat = Atot(ks_valt, Dmin/1e6)
    tcol = t_col(eta, a_pl, M_pl, M_s, areat)
    xpr1 = Xpr1(rho, Dmin, areat, M_s, a_pl, eta, M_pl, L_s)
    xpr2 = Xpr2(areat, M_s, a_pl, eta, M_pl)

    print('dmin ', Dmin)
    print('Dt ', Dt)
    print('Dc ', Dc)

    t_dust1, bmu_dust1, fth_dust1 = F_th_dust(1e-6, L_s, areat, a_pl)
    t_planet1, bmu_planet1, fth_planet1 = F_th_planet(1e-6, L_s, R_pl, a_pl)

    t_dust100, bmu_dust100, fth_dust100 = F_th_dust(1e-4, L_s, areat, a_pl)
    t_planet100, bmu_planet100, fth_planet100 = F_th_planet(1e-4, L_s, R_pl, a_pl)

    Fs_planet1 = F_s(L_s, bmu_planet1, t_planet1, a_pl)
    Fs_scat_planet1 = F_scat(Fs_planet1, R_pl, 0.32, 0.52, a_pl)

    waverange = linspace(1e-7, 0.001, 200)
    T_pl, B_mu_pl, F_th_pl = F_th_planet(waverange, L_s, R_pl, a_pl)
    T_dst, B_mu_dst, F_th_dst = F_th_dust(waverange, L_s, areat, a_pl)

    T_s = stellarTemp(M_s)

    Fs_planet = F_s(L_s, B_mu_pl, T_s, a_pl)
    F_scat_planet = F_scat(Fs_planet, R_pl, 0.32, 0.52, a_pl)

    Fs_dust = F_s(L_s, B_mu_dst, T_s, a_pl)
    F_scat_dust = F_scat(Fs_dust, areat, 0.32, 0.08, a_pl)



    print('D_min = {0:.3e} micro-meters'.format(Dmin))
    print('K_s = {0:.3e}'.format(ks_val))
    print('K_g = {0:.3e}'.format(kg_val))
    print('Area = {0:.3e} au^2'.format(area))
    #print('Rate of collision = {0:.3e} yr^-1'.format(R_cc))
    print('Xc = {0:.3e}'.format(X_c))
    print('Qd = {0:.3e}'.format(Qd))
    print('Hill Radius = {0:.4e} au'.format(rH))
    print('Relative velocity = {0:.4e} ms^-1'.format(vrel))

    print('Rcc = {0:.4e}'.format(R_cc1))
    print('t_nleft = {0:.4e}'.format(tnleft))
    #print('{0:.4e}'.format(((a_pl*eta)**4.13*Qd**0.63*1.7e8*6*M_pl**0.24*rho*Dc)/(1.3e7*n*(4/pi)**2.27*Mtot0)))
    print('n = {0:.3e}'.format(n))
    print('for wavelength = 1 micro-meter')
    print()
    print('F_pl = {0:.3e} Jy'.format(Fs_scat_planet1*1e26))
    print('F_dust = {0:.3e} Jy'.format(fth_dust1*1e-26))
    print()
    print('for wavelength = 100 micro-meter')
    print()
    print('F_pl = {0:.3e} Jy'.format(fth_planet100*1e26))
    print('F_dust = {0:.3e} Jy'.format(fth_dust100*1e-26))

    print('effective temp = {0:.3e}'.format(t_planet1))
    print('B_nu = {0:.3e}'.format(bmu_planet1))

    drange = linspace(100, 250000, 2500)
    #nrange = n_kg(kg_val, drange)
    timerange = logspace(0, 10, 20)
    timerangea = logspace(0, 10, 1000)

    kg_range = []
    area_range = []
    Mtot00 = Mtot0
    for i in range(len(timerange)):
        Mtot_tt = Mtott(Mtot00, R_cc1, timerange[i], 1e111, 0, 0)
        #Dctt = Dct(Dmax, 1e6*tnleft, timerange[i])
        kg_t = kg(Mtot00, Dmax, rho)

        kg_range.append(kg_t)
        Mtot00 = Mtot_tt
        #Dc00 = Dctt


    Mtot00 = Mtot0
    #Dc00 = Dc
    Mtot_t = Mtot0
    Rccrange = []


    tnleftrange = []
    for i in range(len(timerangea)):
        if timerangea[i] > tnleft:
            tnleftrange.append(timerangea[i])

    #params = [const, 1e6*tnleft, Dmax]
    y0 = Mtott(Mtot00, R_cc1, tnleft*1e6, 1e111, 0, 0)
    A = y0/Dmax**3
    Dctt = Dct(Dmax, 1e6*tnleft, 4.5e9)
    Mtot_tt = Mtott(Mtot00, R_cc1, 4.5e9,1e6*tnleft, A, Dctt)
    kg_t = kg(Mtot_tt*5.97e24, Dctt, rho)
    ks_tk = ks1SI(Dt, kg_t)
    #ks_t = ks1(Dt, kg_t)
    number = N_SI(ks_tk, Dmin/1e6)
    sigtot = Atot(ks_tk, Dmin/1e6)
    dens = number/(4*pi*0.4**2*0.2*0.866*(0.35*1.5e8)**3)

    areattt_range = []

    print("dens = {0:.3e}".format(dens))
    print("number = {0:.3e}".format(number))
    print("sigtot = {0:.3e}".format(sigtot/1.5e11**2))
    print("A = {0:.3e}".format(A))
    #soln = integrate.odeint(g, y0, tnleftrange, args=(params,))

    for i in range(len(timerangea)):
        Dctt = Dct(Dmax, 1e6*tnleft, timerangea[i])
        Mtot_tt = Mtott(Mtot00, R_cc1, timerangea[i],1e6*tnleft, A, Dctt)
        Mtot_ttt = Mtott(Mtot00, R_cc1, timerangea[i], 1e111, 0, 0)

        kg_tt = kg(Mtot_ttt, Dmax, rho)
        ks_tt = ks1(Dt, kg_tt)
        areattt = Atot(ks_tt, Dmin/1e6)

        areattt_range.append(areattt)

        kg_t = kg(Mtot_tt, Dctt, rho)
        ks_t = ks1(Dt, kg_t)
        areatt = Atot(ks_t, Dmin/1e6)
        #area = Atot(ks_val, (Dmin/1e6))

        Mtot_t = Mtot_tt

        # print()
        # print(timerangea[i])
        # print()
        # print('Dctt = {0:.3e}'.format(Dctt))
        # print('Dmax = {0:.3e}'.format(Dmax))
        # print('Mtot0 = {0:.3e}'.format(Mtot00))
        # print('Mtot_tt = {0:.3e}'.format(Mtot_tt))
        # print('kg_t = {0:.3e}'.format(kg_t))
        # print('ks_t = {0:.3e}'.format(ks_t))
        # print('ks_val = {0:.3e}'.format(ks_val))
        # print('areatt = {0:.3e}'.format(areatt))
        area_range.append(areatt)
        #Mtot00 = Mtot_tt
        #Dc00 = Dctt
    print("gate1")
    print(area_range)
    #print(soln)
    #soln = array(soln)
    soln_right = []
    #for i in range(len())
    #soln = soln.transpose(soln)
    area_range = array(area_range)
    #a_soln = concatenate(soln, area_range)
    sup_nrange = []
    for i in range(len(kg_range)):
        sup_nrange.append(n_kg(kg_range[i], drange))

    plt.figure(5)
    plt.title('area')
    plt.loglog(timerangea, area_range)
    plt.loglog(timerangea, areattt_range)
    plt.xlim([2e5, 1e10])
    plt.ylim([5e-9, 1e-5])
    plt.show()

    # plt.figure(4)
    # plt.subplot(211)
    # plt.title('F scattered planet')
    # plt.loglog(waverange, F_scat_planet)
    # plt.subplot(212)
    # plt.title('F thermal planet')
    # plt.loglog(waverange, F_th_pl)
    # plt.tight_layout()
    # plt.show()
    #
    # plt.figure(7)
    # plt.subplot(211)
    # plt.title('F scattered dust')
    # plt.loglog(waverange, F_scat_dust)
    # plt.subplot(212)
    # plt.title('F thermal dust')
    # plt.loglog(waverange, F_th_dst)
    # plt.tight_layout()
    # plt.show()

    plt.figure(2)
    plt.title('number of objects')
    for i in range(len(sup_nrange)):
        plt.loglog(drange, sup_nrange[i], label=repr(i))
        plt.xlim([5000, 200000])
        plt.ylim([1e-5, 1e-2])
        plt.legend()
    plt.show()

    Dcttest = Dct(Dmax, 1e6*tnleft, 4.5e9)
    Mtot_ttest = Mtott(Mtot00, R_cc1, 4.5e9, 1e6*tnleft, A, Dcttest)
    Mtot_tttest = Mtott(Mtot00, R_cc1, 4.5e9, 1e111, 0, 0)

    kg_ttest = kg(Mtot_tttest, Dmax, rho)
    ks_ttest = ks1(Dt, kg_ttest)
    areatttest = Atot(ks_ttest, Dmin/1e6)

    kg_test = kg(Mtot_ttest, Dcttest, rho)
    ks_test = ks1(Dt, kg_test)
    areattest = Atot(ks_test, Dmin/1e6)

    print("t = 4.5e9 years")
    print("Rcc = {0:.3e} yr^-1".format(R_cc1))
    print("Mass without correction = {0:.3e} Earth-mass".format(Mtot_tttest))
    print("Mass with correction = {0:.3e} Earth-mass".format(Mtot_ttest))


if __name__ == '__main__':
    main(0.001*0.012345679, 100, 150000, 1500, 0.4, 1, 317.46, 1, 5.2, 6.9911e7) #Jupiter
    #main(0.01*0.012345679, 0.1, 150, 1500, 0.4, 1, 1, 1, 1, 6.4e7) #earth
    #main(0.01*0.012345679, 0.1, 250, 1500, 0.3, 1, 95.16, 1, 10.1, 5.8232e7) #Saturn
