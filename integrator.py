import numpy as np
import random as rand
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy
import sys
import time

xi = []
for i in range(100):
    xi.append([rand.random(), rand.random()])

print xi
def func(pts):
    x = pts[0]
    y = pts[1]
    m = 100
    nMin = -10
    nMax = 10
    k = 10.

    total = 0.
    for i in range(m):
        for j in range(nMax - nMin + 1):
            for l in range(nMax - nMin + 1):
                val = np.exp(-k * ( (x-xi[i][0]+(nMin+j) )**2 + (y-xi[i][1] + (nMin + l))**2 ))
                total += val
    total = total *  k/np.pi / m
    return total

class Region:
    def __init__(self, parent, loc, size, weight):
        self.parent = parent
        self.loc = loc
        self.size = size
        self.weight = weight
        self.points = []
        self.funcVal = []

    def split(self):
        w1 = self.weight/2
        w2 = w1
        
        largest = [0, 0]
        for i in range(len(self.size)):
            if (self.size[i] > largest[0]):
                largest[0] = self.size[i]
                largest[1] = i

        sizeNew = copy.copy(self.size)
        sizeNew[largest[1]] = self.size[largest[1]]/2

        loc1 = copy.copy(self.loc)
        loc2 = copy.copy(self.loc)
        loc2[largest[1]] += sizeNew[largest[1]]

        d1 = Region(self, loc1, sizeNew, w1)
        d2 = Region(self, loc2, sizeNew, w2)

        return d1, d2

    def generatePoint(self):
        D = len(self.size)

        point = []
        for i in range(D):
            r = self.loc[i] + self.size[i] * rand.random()
            point.append(r)

        self.points.append(point)

    def updateWeight(self, val):
        return Region(self.parent, self.loc, self.size, val)


dim = 2
N = 100
repeats = 100
regions = []
startN = 16

def initialize():
    a = Region(0, [0,0], [1.,1.], 1.)
    b,c = a.split()
    d,e = b.split()
    f,g = c.split()

    first = []
    first.append(d.split())
    first.append(e.split())
    first.append(f.split())
    first.append(g.split())

    second = []
    for i in range(len(first)):
        second.append(first[i][0].split())
        second.append(first[i][1].split())
    
    third = []
    for i in range(len(second)):
        third.append(second[i][0])
        third.append(second[i][1])
    regions.append(third)

def generateBatch(c):
    for n in range(N):
        choose = rand.random()
        at = 0
        choice = 0
        for i in range(len(regions[c])):
            reg = regions[c][i]
            at += reg.weight
            if (choose < at):
                choice = i
                break
        regions[c][choice].generatePoint()

def updateWeights(curr):
    w = []
    for i in range(len(regions[curr])):
        wi = 0
        pts = regions[curr][i].points
    
        vol = 1
        for j in range(dim):
            vol *= regions[curr][i].size[j]

        for j in range(len(regions[curr][i].points)):
            funcVal = func(pts[j])
            wi += funcVal**2
            regions[curr][i].funcVal.append(funcVal)

        if (len(pts) > 0):
            wi = (wi * vol**2 / len(pts)) ** (0.5)
        else:
            wi = regions[curr][i].weight
        w.append(wi)

    w_tot = sum(w)
    for i in range(len(w)):
        w[i] = w[i]/w_tot

    highest = updateRegions(curr, w)     
    return highest  

def updateRegions(c, w):
    highest = [0.,0] 
    for i in range(len(w)):
        if (w[i] > highest[0]):
            highest[0] = w[i]
            highest[1] = i
    toAdd = []
    for i in range(len(w)):
        regs = regions[c]
        if (i != highest[1]):
            newReg = Region(regs[i].parent, regs[i].loc, regs[i].size, w[i])
            toAdd.append(newReg)
        else:
            preReg = Region(regs[i].parent, regs[i].loc, regs[i].size, w[i])
            n1, n2 = preReg.split()
            toAdd.append(n1)
            toAdd.append(n2)
    
    regions.append(toAdd)
    return highest[0]

def calcEff(c, h):
    m = len(regions[c]) + startN
    eff = 1./(m*h)
    return eff
    
def calcIntegral(c):
    n = N*c
    M = 0
    for i in range(len(regions)-1):
        M += len(regions[i])-startN + 1
      
    S1 = 0.
    S2 = 0.
    S3 = 0.
    S4 = 0.
    for i in range(len(regions)):
        Step1 = 0.
        Step2 = 0.
        Step3 = 0.
        Step4 = 0.
        for j in range(len(regions[i])):
            r = regions[i][j]

            for k in range(len(r.points)):
                vol = 1
                for l in range(dim):
                    vol *= r.size[l]

                funcVal = r.funcVal[k]
                Step1 += (funcVal * vol / r.weight)
                Step2 += (funcVal * vol / r.weight) ** 2
                Step3 += (funcVal * vol / r.weight) ** 3
                Step4 += (funcVal * vol / r.weight) ** 4

        S1 += Step1 * (i+1)/M *c
        S2 += Step2 * ((i+1)/M) ** 2 * c
        S3 += Step3 * ((i+1)/M) ** 3 * c
        S4 += Step4 * ((i+1)/M) ** 4 * c

    E1 = S1/n
    E2 = (S2 - S1**2./(n))/(n*(n-1.))
    E4 = ((n*(n-1.))/(n*n * n*(n-1.)*(n-2.)*(n-3.)) * (S4 - 4.*S3*S2/n + 3.*S2**2./n) + (1./(n*(n-1.))**2 - 1./(n*(n-1.)*(n-2.)*(n-3.)))*((S2 - S1**2./n)**2.))


    print E1, "+-", (-E2)**(0.5), "+-", (-E4)**(0.25)
    return E1, E2, E4


def plot(c):
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    
    for i in range(len(regions[c])):
        r = regions[c][i]
        ax.add_patch(patches.Rectangle(r.loc, r.size[0], r.size[1], fill=False))
    plt.title("Visualization of integration regions")
    fig.savefig('test.png')

    plt.show()

def run():
    curr = 0

    last_eff = 0
    for i in range(repeats):
        generateBatch(curr)
        highest = updateWeights(curr)
        curr += 1
    E1, E2, E4 = calcIntegral(curr)
    plot(curr)


start_time = time.time()
initialize()
run()
print "Time: ", time.time() - start_time

"""
x = []
y = []
z = []
for i in range(80):
    for j in range(80):
        xNew = float(i)/80
        yNew = float(j)/80
        x.append(xNew)
        y.append(yNew)
        z.append(func([xNew, yNew]))

plt.hexbin(x,y, C=z, gridsize=50)
cb = plt.colorbar()
plt.savefig('pls.png')
"""
