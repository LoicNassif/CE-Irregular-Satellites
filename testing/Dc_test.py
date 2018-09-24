# Test behaviour of equation 20 from Kennedy

from numpy import logspace
import matplotlib.pyplot as plt

def Dc_equation(Dmax, tnleft, alpha, t):
    denominator = (1 + 0.4 * (t - tnleft)/tnleft)**alpha
    return Dmax / denominator

def main():
    t = logspace(0, 10, 25)
    Dmax = 1.5e5
    alpha = 1./1.2
    tnleft = [0.5e8, 1.0e8, 1.5e8, 2.0e8, 2.5e8]
    colour = ['r', 'b', 'y', 'g', 'c']
    for j in range(len(tnleft)):
        Dc_list = []
        for i in range(len(t)):
            Dc_list.append(Dc_equation(Dmax, tnleft[j], alpha, t[i]))

        #make string to display
        #s = "Dc values in meters \n"
        #for i in range(len(t)):
            #s += "{0:.5e}".format(Dc_list[i]) + "\n"

        plt.loglog(t, Dc_list, colour[j]+'o', label="tnleft = {0:.3e}".format(tnleft[j]))
        #plt.text(1e1, 1e4, s, fontsize=8)
        plt.legend()
        #plt.title("At tnleft {0:.3e}".format(tnleft))
        plt.xlabel("time [yr]")
        plt.ylabel("Dc [m]")

    plt.show()
if __name__ == '__main__':
    main()
