#! /usr/bin/env python

'''
A modification of the truss analyser. This one performs simulated annealing to
attempt to improve the bridge design (purely by moving joint positions, not by
changing any member connections)

Copyright (c) 2013 Bennett Lovelady

pseudocode:

s = s0, e = E(s)
s_best = s, e_best = e
k = 0
while k<k_max && e>e_max:    # i.e. time left && not good enough
    T = Temperature(k/k_max)
    s_new = SomeNeighbour(s)
    e_new = E(s_new)
    if P(e,e_new,T) > rand():
        s = s_new, e = e_new
    if e_new < e_best:
        s_best = s_new, e_best = e_new
    k++
return s_best

E(state) returns the energy of a given state
Temperature(float) returns a float to represent the cooling temperature
SomeNeighbour(state) returns a nearby state (e.g. jiggled a bit)
P(e,e',T) returns probability of moving from s->s' at temp T
    usually defined as 1 if e'<e, or exp(-(e'-e)/T) otherwise
Can define restart(): s=s_best, e=e_best to help with cooling

In a truss context the energy is replaced by maximum load which
is maximised rather than minimised

new cooling schedule:
iterate(k++)
if k = kmax: k0 = (kmax+k0)/2; k = k0    # heat it to half of the last "hot" temp
if no improvements in last "m" tries: reset to the best
if "m" resets in a row: convergence achieved! 
'''

import numpy as np
import copy

bridgenumber = raw_input('input the numbered file in "/anneal" containing the .in files: ')
inpath = '/home/bennett/Documents/ensc1002/anneal/' + bridgenumber + '/'
collating = True
outpath = inpath + 'results/'
outputname = str(bridgenumber)
html_output = True


'''
joints: 
[joint #] [x coord (mm)] [y coord (mm)] [degree of freedom]

memers:
[member #] [joint i] [joint j] [dx] [dy] [L] [lij] [mij] [maxL_y] [maxL_b] [maxL]

supports:
[joint # with x-reaction]
[joint # with y-reaction]
[joint # with y-reaction]

applied:
[joint #] [F_x] [F_y]

material:
list of material properties accessed from a dictionary
'''

# ------ a single data structure for the truss ------
class Truss:
    def __init__(self, joints_, members_, supports_, applied_, material_):
        self.joints = joints_ 
        self.members = members_
        self.supports = supports_
        self.applied = applied_
        self.material = material_
        self.C = None
        self.P = None
        self.Q = None
        self.maxL = None
        self.maxL_m = None
        self.singular = None

    def copy(self):
        return copy.deepcopy(self)

    def LoadData(self):
        # applied forces
        n = self.joints.shape[0]
        self.P = np.zeros(2*n)
        num = self.applied.shape[0]-1
        for i in range(1,num+1):
            self.P[2*self.applied[i,0]-2] = -self.applied[i,1]
            self.P[2*self.applied[i,0]-1] = -self.applied[i,2]
        # members table
        for k in range(self.members.shape[0]):
            i   = self.members[k,1]
            j   = self.members[k,2]
            dx  = self.joints[j-1,1] - self.joints[i-1,1]
            dy  = self.joints[j-1,2] - self.joints[i-1,2]
            L   = np.sqrt(dx**2 + dy**2)
            lij = dx/L
            mij = dy/L
            self.members[k,3], self.members[k,4] = dx,dy
            self.members[k,5], self.members[k,6], self.members[k,7] = L, lij, mij
        # coefficient matrix
        Cm = np.zeros((2*n,2*n))
        for k in range(self.members.shape[0]):
            i   = self.members[k,1]
            j   = self.members[k,2]
            lij = self.members[k,6]
            mij = self.members[k,7]    
            Cm[2*i-2, k] = lij
            Cm[2*i-1, k] = mij
            Cm[2*j-2, k] = -lij
            Cm[2*j-1, k] = -mij
        # now add the coefficients for the reaction forces:
        Cm[2*self.supports[0]-2, 2*n-3] = 1
        Cm[2*self.supports[1]-1, 2*n-2] = 1
        Cm[2*self.supports[2]-1, 2*n-1] = 1
        self.C = Cm
        # degrees of freedom for each joint
        # roadbed can move l/r
        for i in range(self.joints.shape[0]):
            if self.joints[i,2] == 0:
                self.joints[i,3] = 1
            else:
                self.joints[i,3] = 3
        # supports can't move
        self.joints[self.supports[0]-1,3] = 0
        self.joints[self.supports[1]-1,3] = 0
        self.joints[self.supports[2]-1,3] = 0
        # loaded joints can't move
        for i in range(1,num+1):
            self.joints[self.applied[i,0]-1,3] = 0

    def Solve(self):
        self.singular = False
        try:
            # invert the matrix
            self.Q = np.linalg.solve(self.C,self.P)
        except np.linalg.LinAlgError:
            # if inversion doesn't work use numerical approximation
            print "singular matrix. using least squares.."
            self.singular = True
            self.Q = np.linalg.lstsq(self.C,self.P)[0]
        for i in range(len(self.Q)):
            # it sometimes gives tiny values like 1e-13, let's get rid of those:
            if abs(self.Q[i]) < 0.01:
                self.Q[i] = 0.0
        # find the maximum load
        for k in range(self.members.shape[0]):
            if self.Q[k] != 0:
                # yield: F/A <= s.f * s_y
                L_yield = self.material['sf']*self.material['s_y']*self.material['A']/abs(self.Q[k])
                self.members[k,8] = L_yield
                # buckling: F <= s.f * pi^2 * E * I / l^2
                if self.Q[k] < 0:
                    L_buckling = self.material['sf']*np.pi**2*self.material['E']*self.material['I'] \
                                    / (abs(self.Q[k])*self.members[k,5]**2)
                    self.members[k,9] = L_buckling
                    # which is the smaller limit?
                    self.members[k,10] = min(self.members[k,8], self.members[k,9])
                else:
                    self.members[k,9] = 0
                    self.members[k,10] = self.members[k,8]
        self.maxL = 1e6
        for k in range(self.members.shape[0]):
            m = self.members[k,10]
            if m != 0 and m <= self.maxL:
                self.maxL = m
        self.maxL_m = []
        for k in range(self.members.shape[0]):
            # to account for floating point errors,
            # check if the member's maxL is within 1% of maxL
            if np.allclose([self.members[k,10]], [self.maxL], rtol=1e-2):
                self.maxL_m.append(k)

    # maximum load of each member as strings
    def MaxLStrings(self):
        maxL_str = ['']*self.members.shape[0]
        for k in range(self.members.shape[0]):
            m = self.members[k,10]
            maxL_str[k] = '---' if m == 0 else str(round(m,2))
        return maxL_str
    
    # maximum x-coord => bridge length
    def TrussLength(self):
        maxX = 0
        for i in range(self.joints.shape[0]):
            if self.joints[i,1] > maxX:
                maxX = self.joints[i,1]
        return maxX
    
    # total length of material used
    def TotalMaterial(self):
        sumL = 0
        for k in range(self.members.shape[0]):
            sumL = sumL + self.members[k,5]
        return sumL
        


# ------ functions for annealing ------


def Temperature(k, kmax):
    t1 = 1-k/kmax
    t2 = np.exp(-k/kmax)
    t3 = 1/(1+np.exp(k/kmax))
    t4 = 2/(1+np.exp(10*k/kmax))
    t5 = 2/(1+np.exp(15*k/kmax))
    return t5

def Prob(load1, load2, temp):
    if load2 > load1:
        return 1.0
    else:
        return np.exp(-10*(load1-load2)/temp)


def Jiggle(jin):
    # takes a "joints" structure and jiggles them a bit
    # dof is 0,1,2,3 = none, l/r, u/d, u/d/l/r movement allowed
    mag = 1.0
    jout = np.zeros((jin.shape[0],4))
    for i in range(jin.shape[0]):
        for j in range(4):
            jout[i,j] = jin[i,j]
    i = int(np.random.random()*jin.shape[0])
    dof = jin[i,3]
    dx, dy = np.random.random(), np.random.random()
    dx, dy = dx*2-1, dy*2-1
    if dof == 0:
        jout[i,1], jout[i,2] = jin[i,1], jin[i,2]
    elif dof == 1:
        jout[i,1], jout[i,2] = jin[i,1] + dx*mag, jin[i,2]
    elif dof == 2:
        jout[i,1], jout[i,2] = jin[i,1], jin[i,2] + dy*mag
    else:
        jout[i,1], jout[i,2] = jin[i,1] + dx*mag, jin[i,2] + dy*mag
    return jout
        

# ------ functions for output ------
# some things used by both regular and html output sections:
# tension as a string
def TensionState(x):
    if (x<0):
        return 'comp.'
    elif (x>0):
        return 'tens.'
    else:
        return ''


def OutputPlaintext():
    maxL_str = MaxLStrings()
    fout = open(outpath + outputname + '.out', 'w')
    if singular:
        fout.write("The system produced a singular matrix.\n\
	            These results were found by least squares approximation.\n")
    # output a table of all members' internal forces
    fout.write('mem. #\tstate\tforce/L\tmax load/N\n')
    for k in range(members.shape[0]):
        fout.write(str(k+1) + '\t' + TensionState(Q[k]) + \
                   '\t' + str(Q[k]) + '\t' + maxL_str[k] + '\n')
    for i in range(3):
        fout.write(reactions[i] + '\t\t' + str(Q[members.shape[0]+i]) + '\n')
    fout.write('\n')
    fout.write('max load:\t' + str(maxL) + ' N\n')
    fout.write('\t\t' + str(maxL/9.81) + ' kg\n\n')
    fout.write('weakest members:\t')
    for k in range(len(maxL_m)):
        fout.write(str(maxL_m[k]+1))
        if k != len(maxL_m)-1:
            fout.write(', ')
        else:
            fout.write('\n\n')
    # now the bridge length and total length of all members
    fout.write('bridge length:\t' + str(int(TrussLength())) + 'mm\n')
    fout.write('total material:\t' + str(int(TotalMaterial())) + 'mm\n')
    fout.close()


def initTable(fout,n):
    fout.write('<table>\n\t<colgroup>\n' + \
                n*'\t\t<col span="1" width="100">\n' + \
               '\t</colgroup>\n')


def OutputHTML(truss_, filenumber):
    print "outputting truss with maxL = " + str(truss_.maxL)
    maxL_str = truss_.MaxLStrings()
    fout = open(outpath + outputname + '.' + str(int(filenumber)) + '.html', 'w')
    fout.write('<!DOCTYPE html>\n<html>\n<head>\n<title>Results</title>\n</head>\n')
    fout.write('\n<body>\n')
    fout.write('<img src="' + outputname + '.' + str(int(filenumber)) + '.png">\n')
    if truss_.singular:
        fout.write('<p>The system produced a singular matrix.<br>\n' + \
	           'These results were found by least squares approximation.<br></p>\n\n')

    # table of internal forces
    fout.write('<p>Internal forces:</p>\n')
    initTable(fout,5)
    fout.write('\t<tr><td>Member #</td><td>Length/mm</td><td>State</td>' + \
               '<td>Force/Load</td><td>Max Load/N</td></tr>\n')
    for k in range(truss_.members.shape[0]):
        fout.write('\t<tr><td>' + str(k+1) + \
                   '</td><td>' + str(round(truss_.members[k,5],1)) + \
                   '</td><td>' + TensionState(truss_.Q[k]) + \
                   '</td><td>' + str(round(truss_.Q[k],4)) + \
                   '</td><td>' + maxL_str[k] + '</td></tr>\n')
    for i in range(3):
        fout.write('\t<tr><td>' + reactions[i] + \
                   '</td><td>' + '---' + \
                   '</td><td>' + str(round(truss_.Q[i-3],4)) + '</td><td></td></tr>\n')
    fout.write('</table>\n\n<br><br>')

    # weakest members
    fout.write('<p>Weakest members: </p>')
    initTable(fout,3)
    fout.write('\t<tr><td>Member</td><td></td><td>Max Load / N</td></tr>\n')
    for i in range(len(truss_.maxL_m)):
        fout.write('\t<tr><td>Member ' + str(truss_.maxL_m[i]+1) + ':' + \
                   '</td><td>' + \
                   '</td><td>' + maxL_str[truss_.maxL_m[i]] + '</td></tr>\n')
    fout.write('</table><br><br>\n\n')

    # maximum load, bridge length and material used
    initTable(fout,3)
    fout.write('\t<tr><td>Max load:</td><td></td><td>' + \
                str(round(truss_.maxL,4)) + ' N</td></tr>\n')
    fout.write('\t<tr><td></td><td></td><td>' + \
                str(round(truss_.maxL/9.81,4)) + ' kg</td></tr>\n')
    fout.write('\t<tr><td> </td></tr>\n')
    fout.write('\t<tr><td>Bridge length:</td><td></td><td>' + \
                str(int(truss_.TrussLength())) + 'mm</td></tr>\n')
    fout.write('\t<tr><td>Material used:</td><td></td><td>' + \
                str(int(truss_.TotalMaterial())) + 'mm</td></tr>\n')
    fout.write('</table><br><br>\n\n')
    fout.close()


def Draw(truss_, filenumber):
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import matplotlib.text as text
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    maxX = truss_.TrussLength()
    unit = maxX/8

    # draw the members first
    for k in range(truss_.members.shape[0]):
        i, j = truss_.members[k,1], truss_.members[k,2]
        ix, iy = truss_.joints[i-1,1], truss_.joints[i-1,2] 
        jx, jy = truss_.joints[j-1,1], truss_.joints[j-1,2]
        c = "#ff0000" if truss_.Q[k]>0 else ("#0000ff" if truss_.Q[k]<0 else "#00ff00")
        ax.add_line(mlines.Line2D([ix,jx], [iy,jy], lw=2., color=c))
        mx = (ix+jx)/2
        my = (iy+jy)/2
        ax.text(mx,my,str(int(truss_.members[k,0])), fontsize=8)

    # now the joints
    for i in range(truss_.joints.shape[0]):
        ax.plot(truss_.joints[i,1], truss_.joints[i,2], 'ko')

    # draw the reaction and load forces. this gets a bit hairy..
    r = truss_.Q.tolist()[-3:]
    for i in range(3):
        xp, yp = truss_.joints[truss_.supports[i]-1,1], truss_.joints[truss_.supports[i]-1,2]
        axis = [1,0] if i == 0 else [0,1]
        axis = [axis[0]*unit*1.2, axis[1]*unit]
        ax.arrow(xp-axis[0],yp-axis[1],axis[0]*0.6,axis[1]*0.6,\
                 fc='k', ec='k', head_width=unit*0.15, head_length=unit*0.25)
        ax.text(xp-axis[0]+unit*0.1, yp-axis[1]+unit*0.1, reactions[i], fontsize=8)
    num = truss_.applied.shape[0]-1
    for i in range(1,num+1):
        xp, yp = truss_.joints[truss_.applied[i,0]-1,1], truss_.joints[truss_.applied[i,0]-1,2]
        lx, ly = truss_.applied[i,1], truss_.applied[i,2]
        mag = np.sqrt(lx**2 + ly**2)
        lx, ly = lx/mag, ly/mag
        ax.arrow(xp+lx*unit*0.2, yp+ly*unit*0.2, lx*unit*0.75, ly*unit*0.75,\
                 fc='k', ec='k', head_width=unit*0.15, head_length=unit*0.25)
        ax.text(xp+(lx*1.5)*unit, yp+(ly*1.5)*unit, 'Load', fontsize=8)
    # draw a scale in the lower left
    xp = (truss_.joints[truss_.supports[0]-1,1] + truss_.joints[truss_.applied[1,0]-1,1])/2 - unit
    yp = truss_.joints[truss_.supports[0]-1,2] - unit*1.5
    ax.add_line(mlines.Line2D([xp,xp+unit],[yp,yp], lw=1., color="k"))
    ax.add_line(mlines.Line2D([xp,xp],[yp-0.1*unit,yp+0.1*unit], lw=1., color="k"))
    ax.add_line(mlines.Line2D([xp+unit,xp+unit],[yp-0.1*unit,yp+0.1*unit], lw=1., color="k"))
    ax.text(xp+unit*0.1, yp+unit*0.1, str(int(unit))+'mm', fontsize=8)
    # write the max load below the scale
    ax.text(xp, yp-0.3*unit, 'max:'+str(round(truss_.maxL,2))+' N', fontsize=8)
    # make sure the truss appears in the centre of the image
    ax.axis('equal')
    ax.axis([-maxX/10.0,maxX*1.1,-maxX/1.5,+maxX/1.5])
    ax.axis('off')
    fig.savefig(outpath + outputname + '.' + str(int(filenumber)) + '.png')



# ------ load initial data ------

# material properties
material = {}
with open(inpath+'material.in', 'r') as fin:
    for line in fin:
        if (line[0] != '#'):
            spl = line.split()
            if len(spl)==2:
                (key, value) = spl 
                material[key] = float(value)

# joints
with open(inpath+'joints.in','r') as fin:
    joints = np.genfromtxt(fin, comments="#", delimiter="\t")

# supports 
with open(inpath+'supports.in', 'r') as fin:
    supports = np.genfromtxt(fin, comments="#")
reactions = ['R'+str(int(supports[0]))+'x',\
             'R'+str(int(supports[1]))+'y',\
             'R'+str(int(supports[2]))+'y']
coord = ['x', 'y', 'y']

# applied forces
with open(inpath+'applied.in','r') as fin:
    applied = np.genfromtxt(fin, comments="#", delimiter="\t")
n = joints.shape[0]
P = np.zeros(2*n)
num = applied.shape[0]-1
for i in range(1,num+1):
    P[2*applied[i,0]-2] = -applied[i,1]
    P[2*applied[i,0]-1] = -applied[i,2]

# members
with open(inpath+'members.in', 'r') as fin:
    members = np.genfromtxt(fin, comments="#", delimiter="\t")
for n in range(8):
    members = np.concatenate((members, [[0]]*members.shape[0]), 1)
# members table now has five extra columns (for geometric data) 
# and three more for max load results

truss = Truss(joints, members, supports, applied, material)
truss.LoadData()
truss.Solve()
OutputHTML(truss, 0)
Draw(truss, 0)


# ------ the main annealing loop ------
# always chooses a better bridge, may choose a worse bridge
# if the temp is high enough
# if it reaches "maxResets" resets in a row, the solution has converged
k0 = 0.0
k = k0
kmax = float(raw_input("length of heating cycle? "))
iterations = 0
fails = 0
failsPerReset = 400
resets = 0
maxResets = 5
truss_best = truss
hitsSinceLastDraw = 0
hitsPerDraw = 20
f = open(outpath + outputname + '.maxL.out','w')
# draw each improvement?
drawingSteps = True

while resets < maxResets:
    joints_new = Jiggle(truss.joints)
    truss_new = Truss(joints_new, members, supports, applied, material)
    truss_new.LoadData()
    truss_new.Solve()
    if truss_new.maxL > truss_best.maxL:
        truss_best = truss_new.copy()
        resets = 0
        # draw it?
        if drawingSteps:
            hitsSinceLastDraw = hitsSinceLastDraw + 1
            if hitsSinceLastDraw == hitsPerDraw:
                hitsSinceLastDraw = 0
                Draw(truss_best, iterations)
    T = Temperature(k, kmax)
    Pr = Prob(truss.maxL, truss_new.maxL, T)
    if Pr > np.random.random()**2:
        # good stuff, continue
        f.write('# ' + str(int(iterations)) + ':\t' + str(truss_new.maxL) + '\t' \
            + str(T) + '\t' + str(Pr) + '\n')
        truss = truss_new.copy()
    else:
        fails = fails + 1

    iterations = iterations + 1
    if (iterations % 500) == 0:
        print "iteration #" + str(int(iterations)) + " - T: " + str(T)
    
    k = k + 1
    if k >= kmax:
        #k0 = (kmax + k0)/2
        k = k0
        if abs(k-kmax) < 100:
            # don't bother heating by less than 50 units
            print "can't heat anymore.. ending simulation"
            break
        print "heating! k: " + str(k) + ", T: " + str(Temperature(k, kmax))

    if fails == failsPerReset:
        fails = 0
        resets = resets + 1
        truss = truss_best.copy()
        print "..reset! " + str(maxResets-resets) + " resets remain"

f.close()


# ------ write the output with the best configuration ------
print "last truss is maxL = " + str(truss.maxL)
print "best truss is maxL = " + str(truss_best.maxL)
OutputHTML(truss_best, iterations)
Draw(truss_best, iterations)


print "analysis complete!"
