from swarms import CollSwarm

#main()
def main():
    jupiter = CollSwarm(7.37307e19, 100, 150000, 1, 1, 317.46, 5.2, 6.9911e7, 0.4)
    area = jupiter.computeAtot()
    print('area = {0:.3e}'.format(area))

if __name__ == '__main__':
    main()
