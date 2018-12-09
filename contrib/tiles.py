import numpy as np
from matplotlib import pyplot as pl
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree as kdtree
from timeit import timeit

def points(n):
    # getting some random points
    n1 = n
    n2 = 2 * n1
    n3 = 2 * n2
    x = np.random.random(n1)
    y = np.random.random(n1)
    x = np.append(x, np.random.random(n2) * 0.6666)
    y = np.append(y, np.random.random(n2) * 0.6666)
    x = np.append(x, np.random.random(n3) * 0.3333)
    y = np.append(y, np.random.random(n3) * 0.3333)
    return (x, y)

def density(x, y, r = 0.1):
    # getting the point density for each point
    disk = np.pi * r**2
    tree = kdtree(np.transpose((x, y)))
    lsts = tree.query_ball_tree(tree, r = r)
    m = len(lsts)
    dens = np.zeros(m)
    for i in range(m):
        dens[i] = len(lsts[i])
    return dens / disk

def coarser(x, y, r = 0.1):
    # getting the point density for a coarser subset (seeds)
    # of points that are close to grid cell centers (xc, yc)
    cc = np.arange(r/4, 1, r/2)
    xc, yc = np.meshgrid(cc, cc)
    xc, yc = xc.flatten(), yc.flatten()

    # kdtrees for the points and the cell centers
    grid = kdtree(np.transpose((xc, yc)))
    tree = kdtree(np.transpose((x, y)))

    # those points in tree that are close to grid
    seed = grid.query_ball_tree(tree, r = r/2)
    seeds = [x for x in seed if x != []]

    # as in density
    lsts = tree.query_ball_tree(tree, r = r)
    m = len(lsts)
    dens = np.zeros(m)

    # but now a loop only over seeds
    for l in seeds:
        i = l[0]
        dens[i] = len(lsts[i])
    disk = np.pi * r**2
    return (xc, yc, dens / disk)

def show(r = 0.1):
    x, y = points(100)
    #p = density(x, y, r)
    xc, yc, p = coarser(x, y, r)
    phi = np.arange(0, np.pi*2, 0.05)
    pl.plot(r * np.cos(phi), r * np.sin(phi), 'r')
    pl.scatter(xc, yc, s = 1, c = 'r')
    pl.scatter(x, y, s = 20, c = p, norm = LogNorm())
    #pl.scatter(x, y, s = 20, c = p)
    pl.colorbar()
    pl.scatter(x, y, s = 1, c = 'k')
    pl.grid()
    pl.axes().set_aspect('equal')
    pl.show()

def test_density():
    density(x, y)

def test_coarser():
    coarser(x, y)

if __name__ == "__main__":
    show()
    num = 100
    print('timing each %i times' % num)
    x, y = points(100)
    print(timeit('test_density()', setup = 'from __main__ import test_density', number = num))
    print(timeit('test_coarser()', setup = 'from __main__ import test_coarser', number = num))
