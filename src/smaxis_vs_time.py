# Loic
import swarms
from numpy import pi, sqrt, linspace, array
import matplotlib.pyplot as plt

RHO=1500

def a_pl_comp(s, t):
    Qd = s.computeQd(s.Dc)
    part1 = 1.3e7 * (s.M_init/5.972e24) * (s.M_s/1.989e30)**1.38 * s.f_vrel**2.27
    part2 = Qd**0.63 * s.rho * (s.Dmax/1000) * (s.M_pl/5.972e24)**0.24 * s.eta**4.13
    return (t  * part1 / part2)**(1./4.13)

def Fth(s, t, a_pl):
    Qd = s.computeQd(s.Dc)
    part1Rcc = 1.3e7 * (s.M_init/5.972e24) * (s.M_s/1.989e30)**1.38 * s.f_vrel**2.27
    part2Rcc = Qd**0.63 * s.rho * (s.Dmax/1000) * (s.M_pl/5.972e24)**0.24 * s.eta**4.13 * a_pl**4.13
    Rcc = part1Rcc / part2Rcc
    part1 = (1/3.9e-6) * (1/(s.rho * (s.d_pl/1.496e11)**2)) * (s.Dc/1000)**-0.9 * (s.Dmin*1e6)**-0.7
    part2 = (s.M_init/5.972e24) / (1 + Rcc * t)
    return part1 * part2

#main
def main(swarm_argv, lamb, t):
    M0 = swarm_argv[0]; Dt = swarm_argv[1]
    Dmax = swarm_argv[2]; L_s = swarm_argv[3]
    M_s = swarm_argv[4]; M_pl = swarm_argv[5]
    a_pl = swarm_argv[6]; R_pl = swarm_argv[7]
    eta = swarm_argv[8]; Nstr = swarm_argv[9]
    d_pl = swarm_argv[10]

    semi_major = []
    fth = []

    s = swarms.CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl,
                            a_pl, R_pl, eta, Nstr, d_pl,
                            rho=RHO, fQ=5, f_vrel=4/pi,
                            correction=True, alpha=1.2)

    for i in range(len(t)):
        s.updateSwarm(t[i])

        a_plv = a_pl_comp(s, t[i])
        #T = s.computeT(L_s, a_plv)
        #bnu = s.computeBmu(lamb, T)
        s2 = swarms.CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl,
                                a_plv, R_pl, eta, Nstr, d_pl,
                                rho=RHO, fQ=5, f_vrel=4/pi,
                                correction=True, alpha=1.2)
        F_th = s2.computeFth(array([lamb]), swarm=True)
        fth.append(F_th/1e-26)
        semi_major.append(a_plv)

        #print("time: "+str(t[i])+"\t semi-major: "+str(a_plv)+"\t thermal: "+str(F_th))


    plt.figure(1)
    plt.plot(t, semi_major)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.xlabel("time [yr]")
    plt.ylabel("semi-major axis [au]")
    plt.loglog()
    plt.show()

    plt.figure(2)
    plt.plot(t, fth)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.xlabel("time [yr]")
    plt.ylabel("F_th [Jy]")
    plt.loglog()
    plt.show()

# if __name__ == '__main__':
#     M0 = 10 * 7.34767309e22; Dt = 100.; Dmax = 250000.; L_s = 20 * 3.828e26;
#     M_s = 1.86 * 1.989e30; M_pl = 1.89587112e27; a_pl = 7.48e12
#     R_pl = 6.9911e7; eta = 0.4; Nstr = 6.; d_pl = 3.086e17
#
#     t = linspace(1e6, 10e10, 500)
#     lamb = 1e-4
#
#     argv = [M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr, d_pl]
#     main(argv, lamb, t)
