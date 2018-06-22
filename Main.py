from swarms import CollSwarm
from numpy import linspace, zeros
import matplotlib.pyplot as plt

#main()
def main():
    #CollSwarm(self, M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, Nstr, rho=1500)
    jupiter = CollSwarm(7.37307e19, 100., 150000., 3.828e26, 1.989e30,
                        1.89587112e27, 7.7792e11, 6.9911e7, 0.4, 6.)

    jupiter3 = CollSwarm(7.37307e19, 100., 150000., 3.828e26, 1.989e30,
                        1.89587112e27, 7.7792e11, 6.9911e7, 0.4, 6., correction=False)

    time = linspace(0, 1e10, 50000)
    area = zeros(50000)
    area3 = zeros(50000)

    for i in range(len(time)):
        a = jupiter.computeAtot()*4.4683677582981e-23
        b = jupiter3.computeAtot()*4.4683677582981e-23
        area[i] = a
        area3[i] = b
        #print("tnleft = {0:.3e}".format(jupiter.computetnleft(jupiter.Rcc0)))
        #mass = jupiter.computeMtot(time[i])
        #print("mass = {0:.3e}".format(mass))
        jupiter.updateSwarm(time[i])
        jupiter3.updateSwarm(time[i])
        #N = jupiter.computeNtot()
        if time[i] % 100 == 0:
            print("time = {0:.3e}".format(time[i]))
            print('area = {0:.6e} AU^2 \n'.format(a))
            print('area = {0:.6e} AU^2 \n'.format(b))
    #print('number of particles below Dt = {0:.6e}'.format(N))

    jupiter2 = CollSwarm(7.37307e19, 100., 150000., 3.828e26, 1.989e30,
                        1.89587112e27, 7.7792e11, 6.9911e7, 0.4, 6.)

    print("t = 4.5e9 years")
    print("Ntot = {0:.3e}".format(jupiter2.computeNtot()))
    print("Mass = {0:.3e} Earth-mass".format(jupiter2.computeMtot(4.5e9)/5.972e24))
    print("Rcc = {0:.3e} yr^-1".format(jupiter2.Rcc0))
    jupiter2.updateSwarm(4.5e9)
    print("Area = {0:.3e} AU^2".format(jupiter2.computeAtot()*4.4683677582981e-23))
    print("tnleft = {0:.3e} ".format(jupiter2.tnleft/1e6))

    plt.figure(1)
    plt.loglog(time, area, 'r', label="corrected stranding evolution")
    plt.loglog(time, area3, ls='--', label="static Dc")
    plt.axvline(jupiter2.tnleft, color='g', ls='--', label="first stranded object")
    plt.xlabel("time [yr]")
    plt.ylabel("area [AU^2]")
    plt.title("Area of Circumplanetary Swarm")
    plt.legend()
    plt.xlim([2e5, 1e10])
    plt.ylim([5e-9, 1e-5])
    plt.show()

    #t = 4.5e9
    #jupiter.updateSwarm(t)
    #area = jupiter.computeAtot()
    #N = jupiter.computeNtot(dhigh=100)
    #print("At t = 4.5e9 years")
    #print('area = {0:.6e} AU^2'.format(area*4.4683677582981e-23))
    #print('number of particles below Dt = {0:.6e}'.format(N))

if __name__ == '__main__':
    main()
