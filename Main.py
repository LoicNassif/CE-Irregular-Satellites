from swarms import CollSwarm

#main()
def main():
    #CollSwarm(self, M0, Dt, Dmax, L_s, M_s, M_pl, a_pl, R_pl, eta, rho=1500)
    jupiter = CollSwarm(7.37307e19, 100, 150000, 3.828e26, 1.989e30,
                        1.89587112e27, 7.7792e11, 6.9911e7, 0.4)
    area = jupiter.computeAtot()
    print('area = {0:.3e}'.format(area))

if __name__ == '__main__':
    main()
