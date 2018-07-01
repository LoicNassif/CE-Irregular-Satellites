from swarms import CollSwarm
from numpy import linspace, zeros, logspace
import matplotlib.pyplot as plt

#main()
def main():
    # Saturn config
    #M0 = 7.37307e19; Dt = 100.; Dmax = 250000.; L_s = 3.828e26;
    #M_s = 1.989e30; M_pl = 5.683e26; a_pl = 1.433537e12
    #R_pl = 5.8232e7; eta = 0.3; Nstr = 2.

    # Neptune config
    #M0 = 7.37307e19; Dt = 100.; Dmax = 250000.; L_s = 3.828e26;
    #M_s = 1.989e30; M_pl = 1.0243e26; a_pl = 4.50439e12
    #R_pl = 2.4622e7; eta = 0.2; Nstr = 6.

    # Jupiter config
    M0 = 7.37307e19; Dt = 100.; Dmax = 150000.; L_s = 3.828e26;
    M_s = 1.989e30; M_pl = 1.89587112e27; a_pl = 7.7792e11
    R_pl = 6.9911e7; eta = 0.4; Nstr = 6.

    jupiter = CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr)
    jupiter3 = CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr,
                            correction=False)

    time = linspace(0, 1e10, 50000)
    area = zeros(50000)
    area3 = zeros(50000)

    for i in range(len(time)):
        a = jupiter.computeAtot()*4.4683677582981e-23
        b = jupiter3.computeAtot()*4.4683677582981e-23
        area[i] = a
        area3[i] = b
        jupiter.updateSwarm(time[i])
        jupiter3.updateSwarm(time[i])

    plt.figure(1)
    plt.loglog(time, area, 'r', label="corrected stranding evolution")
    plt.loglog(time, area3, ls='--', label="static Dc")
    plt.axvline(jupiter.tnleft, color='g', ls='--', label="first stranded object")
    plt.xlabel("time [yr]")
    plt.ylabel("area [AU^2]")
    plt.title("Area of Circumplanetary Swarm - Jupiter")
    plt.legend()
    plt.xlim([2e5, 9e9])
    plt.ylim([5e-9, 1e-5])
    plt.show()

    jupiter4 = CollSwarm(M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr,
                        correction=True)

    drange = linspace(0.1, 250000, 5000)
    time_log = logspace(0, 10, 20)
    num_distribution = []
    for i in range(len(time_log)):
        jupiter4.updateSwarm(time_log[i])
        num = []
        for j in range(len(drange)):
            #print(drange[i])
            num.append(jupiter4.computen(drange[j]))
        num_distribution.append(num)

    plt.figure(2)
    for i in range(len(num_distribution)):
        plt.loglog(drange, num_distribution[i])
    plt.xlim([5000, 200000])
    plt.ylim([1e-5, 1e-2])
    plt.show()

if __name__ == '__main__':
    main()
