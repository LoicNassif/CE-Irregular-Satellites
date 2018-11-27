# Loic
from numpy import logspace, array, linspace, pi
import swarms
import matplotlib.pyplot as plt

RHO = 1500

def Fstar(Ls, Bnu, Ts, dpl):
    sig = 5.670367e-8 #Stefan-Boltzmann constant
    part1 = Ls * Bnu
    part2 = 4 * sig * Ts ** 4 * dpl ** 2
    return part1 / part2

def main(swarm_argv, lamb, t, type_star, planet=False, Ms=None, Mtot0=None, d_plv=None, a_plv=None):
    M0 = swarm_argv[0]; Dt = swarm_argv[1]
    Dmax = swarm_argv[2]; L_s = swarm_argv[3]
    M_s = swarm_argv[4]; M_pl = swarm_argv[5]
    a_pl = swarm_argv[6]; R_pl = swarm_argv[7]
    eta = swarm_argv[8]; Nstr = swarm_argv[9]
    d_pl = swarm_argv[10]

    if type_star == "A":
        M_s = 2.1 * 1.989e30
        L_s = 20 * 3.828e26
    elif type_star == "G":
        M_s = 1 * 1.989e30
        L_s = 1 * 3.828e26
    elif type_star == "M":
        M_s = 0.21 * 1.989e30
        L_s = 0.0079 * 3.828e26

    s = swarms.CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr,
                        d_pl, correction=True, alpha=1./1.2)
    s.updateSwarm(t)
    T_star = s.stellarTemp()
    B_nu = s.computeBmu(lamb, T_star)
    F_star = Fstar(L_s, B_nu, T_star, d_pl)/1e-26
    print("Fstar")
    print(F_star)

    waverange = logspace(-7, -3, 200)
    B_nu_list = s.computeBmu(waverange, T_star)
    F_star_list = Fstar(L_s, B_nu_list, T_star, d_pl)/1e-26

    fth_list = []
    if (a_plv is not None):
        for i in range(len(a_plv)):
            s2 = swarms.CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl,
                                    a_plv[i], R_pl, eta, Nstr, d_pl,
                                    rho=RHO, fQ=5, f_vrel=4/pi,
                                    correction=True, alpha=1.2)

            s2.updateSwarm(t)

            #T = swarm.computeT(L_s, a_plv[i])
            #bnu = swarm.computeBmu(lamb, T)
            if planet:
                F_th = s2.computeFth(array([lamb]), planet=True)
            else:
                F_th = s2.computeFth(array([lamb]), swarm=True)

            fth_list.append(F_th[0]/1e-26)

    contrast_ratio = array(fth_list) / F_star
    #s.updateSwarm(t)
    if planet:
        Fth_swarm = s.computeFth(waverange, planet=True)/1e-26
        Fs_swarm = s.computeFs(waverange, 0.32, 0.08, planet=True)/1e-26
    else:
        Fth_swarm = s.computeFth(waverange, swarm=True)/1e-26
        Fs_swarm = s.computeFs(waverange, 0.32, 0.08, swarm=True)/1e-26

    # print("Fth")
    # print(Fth_swarm)
    # print("Fs")
    # print(Fs_swarm)

    #contrast_ratio_s = []
    # for i in range(len(Fs_swarm)):
    #     if 10e-11 <= Fs_swarm[i] <= 10e-3:
    #         contrast_ratio_s.append(Fs_swarm[i] / F_star[i])
    #         lamb.append(waverange[i])
    contrast_ratio_th = Fth_swarm / F_star_list
    contrast_ratio_s = Fs_swarm / F_star_list

    plt.figure(1)
    plt.loglog(waverange, Fs_swarm, 'r')
    plt.loglog(waverange, Fth_swarm, 'b')
    plt.loglog(waverange, F_star_list, 'g')
    plt.xlabel("wavelength [m]")
    plt.ylabel("F [Jy]")
    plt.ylim([10e-11, 10e5])
    plt.show()

    plt.figure(2)
    #plt.loglog(waverange, contrast_ratio_s, 'r')
    plt.loglog(waverange, contrast_ratio_th, 'b')
    plt.xlabel("wavelength [m]")
    plt.ylabel("F_swarm / F_star")
    plt.ylim([10e-5, 10])
    plt.show()

    plt.figure(3)
    plt.loglog(a_plv/1.496e11, contrast_ratio)
    plt.xlabel("a_pl [au]")
    plt.ylabel("F_swarm / F_star")
    plt.show()


if __name__ == '__main__':
    M0 = 10 * 7.34767309e22; Dt = 100.; Dmax = 250000.; L_s = 10 * 3.828e26;
    M_s = 1.86 * 1.989e30; M_pl = 1 * 1.89587112e27; a_pl = 50 * 1.496e11
    R_pl = 6.9911e7; eta = 0.4; Nstr = 6.; d_pl = 10 * 3.086e16
    argv = [M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr, d_pl]

    t = 1e7
    lamb = 1e-4
    apl = linspace(1 * 1.496e11, 50 * 1.496e11, 500)

    main(argv, lamb, t, "A", a_plv=apl)
