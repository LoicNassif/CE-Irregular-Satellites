"""Computing and plotting thermal emission"""
import swarms
import matplotlib.pyplot as plt
from numpy import pi, linspace

# Define global vars
RHO = 1500

def Fth(B_nu, Dc, Dmin, Mtot0, Rcc0, t, d_pl):
    part1 = (1/3.9e-6) * (1/(RHO * d_pl)) * Dc**-0.9 * Dmin**-0.7
    part2 = Mtot0 / (1 + Rcc0 * t)
    return part1 * part2

def plots(x, y, title, xlabel, ylabel):
    plt.plot(x, y)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

#main
def main(swarm_argv, lamb, t, Ms=None, Mtot0=None, d_plv=None, a_plv=None):
    M0 = swarm_argv[0]; Dt = swarm_argv[1]
    Dmax = swarm_argv[2]; L_s = swarm_argv[3]
    M_s = swarm_argv[4]; M_pl = swarm_argv[5]
    a_pl = swarm_argv[6]; R_pl = swarm_argv[7]
    eta = swarm_argv[8]; Nstr = swarm_argv[9]
    d_pl = swarm_argv[10]

    if a_plv is not None:
        max_i = len(a_plv)
    else:
        a_plv = [a_pl]
        max_i = 1

    fth_list = []
    for i in range(max_i):
        swarm = swarms.CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl,
                                a_plv[i], R_pl, eta, Nstr, d_pl,
                                rho=RHO, fQ=5, f_vrel=4/pi,
                                correction=True, alpha=1./1.2)

        swarm.updateSwarm(t)

        T = swarm.computeT(L_s, a_plv[i])
        bnu = swarm.computeBmu(lamb, T)
        F_th = Fth(bnu, swarm.Dc, swarm.Dmin, M0, swarm.Rcc0, t, d_pl)
        fth_list.append(F_th)

    #plotting
    if len(a_plv) > 1:
        plots(a_plv/1.496e11, fth_list,
            "thermal emission with respect to semi-major axis",
            "a_pl [au]", "F_th [W sr^-1 m^-3]")
        plt.show()

# if __name__ == '__main__':
#     M0 = 10 * 7.34767309e22; Dt = 100.; Dmax = 250000.; L_s = 20 * 3.828e26;
#     M_s = 1.86 * 1.989e30; M_pl = 1.89587112e27; a_pl = 7.48e12
#     R_pl = 6.9911e7; eta = 0.4; Nstr = 6.; d_pl = 3.086e17
#     lamb = 1e-4
#
#     t = 4.5e9
#
#     argv = [M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr, d_pl]
#     aplv = linspace(20 * 1.496e11, 100 * 1.496e11)
#     main(argv, lamb, t, a_plv = aplv)
