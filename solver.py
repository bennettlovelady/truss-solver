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
'''

import numpy as np

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


# ------ some functions for analysis ------

def LoadMembersTable(joints,members_in):
    # find the geometric properties of all the members
    members = members_in
    for k in range(members.shape[0]):
        i   = members[k,1]
        j   = members[k,2]
        dx  = joints[j-1,1] - joints[i-1,1]
        dy  = joints[j-1,2] - joints[i-1,2]
        L   = np.sqrt(dx**2 + dy**2)
        lij = dx/L
        mij = dy/L
        members[k,3], members[k,4], members[k,5], members[k,6], members[k,7]\
            = dx, dy, L, lij, mij
    return members


def ConstructCoefficientMatrix(joints, members, supports):
    n = joints.shape[0] 
    C = np.zeros((2*n,2*n))
    for k in range(members.shape[0]):
        i   = members[k,1]
        j   = members[k,2]
        lij = members[k,6]
        mij = members[k,7]    
        C[2*i-2, k] = lij
        C[2*i-1, k] = mij
        C[2*j-2, k] = -lij
        C[2*j-1, k] = -mij
    # now add the coefficients for the reaction forces:
    C[2*supports[0]-2, 2*n-3] = 1
    C[2*supports[1]-1, 2*n-2] = 1
    C[2*supports[2]-1, 2*n-1] = 1
    return C


def SolveTruss(C, P):
    singular = False
    try:
        # invert the matrix
        Q = np.linalg.solve(C,P)
    except np.linalg.LinAlgError:
        # if inversion doesn't work use numerical approximation
        print "singular matrix. using least squares.."
        singular = True
        Q = np.linalg.lstsq(C,P)[0]
    for i in range(len(Q)):
        # it sometimes gives tiny values like 1e-13, let's get rid of those:
        if abs(Q[i]) < 0.001:
            Q[i] = 0.0
    return singular, Q


def FindMaxLoad():
    # use the material properties in materials.in to find the safe load
    # go through every member and calculate its maximum safe load
    for k in range(members.shape[0]):
        if Q[k] != 0:
            # yield: F/A <= s.f * s_y
            L_yield = material['sf'] * material['s_y'] * material['A'] / abs(Q[k])
            members[k,8] = L_yield
            # buckling: F <= s.f * pi^2 * E * I / l^2
            if Q[k] < 0:
                members[k,9] = material['sf'] * np.pi**2 * material['E'] * material['I'] \
                                / (abs(Q[k]) * members[k,5]**2)
                # which is the smaller limit?
                members[k,10] = min(members[k,8], members[k,9])
            else:
                members[k,10] = members[k,8]
    maxL = 1e6
    for k in range(members.shape[0]):
        m = members[k,10]
        if m != 0 and m <= maxL:
            maxL = m
    maxL_m = []
    for k in range(members.shape[0]):
        # to account for floating point errors,
        # check if the member's maxL is within 1% of maxL
        if np.allclose([members[k,10]], [maxL], rtol=1e-2):
            maxL_m.append(k)
    return maxL, maxL_m


# ------ functions for annealing ------


def Temperature(k, kmax):
    t1 = 1-k/kmax
    t2 = np.exp(-k/kmax)
    t3 = 1/(1+np.exp(k/kmax))
    t4 = 1.5/(1+np.exp(10*k/kmax))
    return t4

def Prob(load1, load2, temp):
    if load2 > load1:
        return 1.0
    else:
        return np.exp(-(load1-load2)/temp)


def Jiggle(jin):
    # takes a "joints" structure and jiggles them a bit
    # dof is 0,1,2,3 = none, l/r, u/d, u/d/l/r movement allowed
    mag = 1.0
    jout = jin
    for i in range(jout.shape[0]):
        dof = jout[i,3]
        dx, dy = np.random.random(), np.random.random()
        dx, dy = dx*2-1, dy*2-1
        if dof in [1,3]:
            jout[i,1] = jout[i,1] + dx*mag
        if dof in [2,3]:
            jout[i,2] = jout[i,2] + dy*mag
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

# maximum load of each member as strings
def MaxLStrings():
    maxL_str = ['']*members.shape[0]
    for k in range(members.shape[0]):
        m = members[k,10]
        maxL_str[k] = '---' if m == 0 else str(round(m,2))
    return maxL_str

# maximum x-coord => bridge length
def TrussLength():
    maxX = 0
    for i in range(joints.shape[0]):
        if joints[i,1] > maxX:
            maxX = joints[i,1]
    return maxX

# total length of material used
def TotalMaterial():
    sumL = 0
    for k in range(members.shape[0]):
        sumL = sumL + members[k,5]
    return sumL

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


def OutputHTML(filenumber):
    maxL_str = MaxLStrings()
    fout = open(outpath + outputname + '.' + str(int(filenumber)) + '.html', 'w')
    fout.write('<!DOCTYPE html>\n<html>\n<head>\n<title>Results</title>\n</head>\n')
    fout.write('\n<body>\n')
    fout.write('<img src="' + outputname + '.' + str(int(filenumber)) + '.png">\n')
    if singular:
        fout.write('<p>The system produced a singular matrix.<br>\n' + \
	           'These results were found by least squares approximation.<br></p>\n\n')

    # table of internal forces
    fout.write('<p>Internal forces:</p>\n')
    initTable(fout,4)
    fout.write('\t<tr><td>Member #</td><td>State</td>' + \
               '<td>Force/Load</td><td>Max Load/N</td></tr>\n')
    for k in range(members.shape[0]):
        fout.write('\t<tr><td>' + str(k+1) + \
                   '</td><td>' + TensionState(Q[k]) + \
                   '</td><td>' + str(round(Q[k],4)) + \
                   '</td><td>' + maxL_str[k] + '</td></tr>\n')
    for i in range(3):
        fout.write('\t<tr><td>' + reactions[i] + \
                   '</td><td>' + '---' + \
                   '</td><td>' + str(round(Q[i-3],4)) + '</td><td></td></tr>\n')
    fout.write('</table>\n\n<br><br>')

    # weakest members
    fout.write('<p>Weakest members: </p>')
    initTable(fout,3)
    fout.write('\t<tr><td>Member</td><td></td><td>Max Load / N</td></tr>\n')
    for i in range(len(maxL_m)):
        fout.write('\t<tr><td>Member ' + str(maxL_m[i]+1) + ':' + \
                   '</td><td>' + \
                   '</td><td>' + maxL_str[maxL_m[i]] + '</td></tr>\n')
    fout.write('</table><br><br>\n\n')

    # maximum load, bridge length and material used
    initTable(fout,3)
    fout.write('\t<tr><td>Max load:</td><td></td><td>' + \
                str(round(maxL,4)) + ' N</td></tr>\n')
    fout.write('\t<tr><td></td><td></td><td>' + \
                str(round(maxL/9.81,4)) + ' kg</td></tr>\n')
    fout.write('\t<tr><td> </td></tr>\n')
    fout.write('\t<tr><td>Bridge length:</td><td></td><td>' + \
                str(int(TrussLength())) + 'mm</td></tr>\n')
    fout.write('\t<tr><td>Material used:</td><td></td><td>' + \
                str(int(TotalMaterial())) + 'mm</td></tr>\n')
    fout.write('</table><br><br>\n\n')
    fout.close()


def Draw(filenumber):
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import matplotlib.text as text
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    # find the unit length to use as padding
    unit = 1000000
    for k in range(members.shape[0]):
        if members[k,5] < unit:
            unit = members[k,5]
    maxX = TrussLength()

    # draw the members first
    for k in range(members.shape[0]):
        i, j = members[k,1], members[k,2]
        ix, iy, jx, jy = joints[i-1,1], joints[i-1,2], joints[j-1,1], joints[j-1,2]
        c = "#ff0000" if Q[k]>0 else ("#0000ff" if Q[k]<0 else "#00ff00")
        ax.add_line(mlines.Line2D([ix,jx], [iy,jy], lw=2., color=c))
        mx = (ix+jx)/2
        my = (iy+jy)/2
        ax.text(mx,my,str(int(members[k,0])), fontsize=8)

    # now the joints
    for i in range(joints.shape[0]):
        ax.plot(joints[i,1], joints[i,2], 'ko')

    # draw the reaction and load forces. this gets a bit hairy..
    r = Q.tolist()[-3:]
    for i in range(3):
        xp, yp = joints[supports[i]-1,1], joints[supports[i]-1,2]
        axis = [1,0] if i == 0 else [0,1]
        axis = [axis[0]*unit*1.2, axis[1]*unit]
        ax.arrow(xp-axis[0],yp-axis[1],axis[0]*0.6,axis[1]*0.6,\
                 fc='k', ec='k', head_width=unit*0.15, head_length=unit*0.25)
        ax.text(xp-axis[0]+unit*0.1, yp-axis[1]+unit*0.1, reactions[i], fontsize=8)
    num = applied.shape[0]-1
    for i in range(1,num+1):
        xp, yp = joints[applied[i,0]-1,1], joints[applied[i,0]-1,2]
        lx, ly = applied[i,1], applied[i,2]
        mag = np.sqrt(lx**2 + ly**2)
        lx, ly = lx/mag, ly/mag
        ax.arrow(xp+lx*unit*0.2, yp+ly*unit*0.2, lx*unit*0.75, ly*unit*0.75,\
                 fc='k', ec='k', head_width=unit*0.15, head_length=unit*0.25)
        ax.text(xp+(lx*1.5)*unit, yp+(ly*1.5)*unit, 'Load', fontsize=8)
    # draw a scale in the lower left
    xp = (joints[supports[0]-1,1] + joints[applied[1,0]-1,1])/2 - unit
    yp = joints[supports[0]-1,2] - unit*1.5
    ax.add_line(mlines.Line2D([xp,xp+unit],[yp,yp], lw=1., color="k"))
    ax.add_line(mlines.Line2D([xp,xp],[yp-0.1*unit,yp+0.1*unit], lw=1., color="k"))
    ax.add_line(mlines.Line2D([xp+unit,xp+unit],[yp-0.1*unit,yp+0.1*unit], lw=1., color="k"))
    ax.text(xp+unit*0.1, yp+unit*0.1, str(int(unit))+'mm', fontsize=8)
    # make sure the truss appears in the centre of the image
    ax.axis('equal')
    ax.axis([-unit,maxX+unit,-maxX/1.5,+maxX/1.5])
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
members = LoadMembersTable(joints, members)


# ------ initial calculation ------
C = ConstructCoefficientMatrix(joints, members, supports)
singular, Q = SolveTruss(C, P) 
maxL, maxL_m = FindMaxLoad()
OutputHTML(0)
Draw(0)


# ------ data is loaded, begin annealing ------
iterations = 0.0
maxiterations = 1000.0
maxL_best = 0.0
f = open(outpath + outputname + '.maxL.out','w')


while iterations < maxiterations:
    joints_new = Jiggle(joints)
    members = LoadMembersTable(joints_new, members)
    C = ConstructCoefficientMatrix(joints_new, members, supports)
    singular, Q = SolveTruss(C, P)
    maxL_new, maxL_m = FindMaxLoad()
    if maxL_new > maxL_best:
        joints_best = joints_new
        maxL_best = maxL_new
    T = Temperature(iterations, maxiterations)
    if Prob(maxL, maxL_new, T) > np.random.random()**2:
        # good stuff, continue
        #Draw(iterations)
        f.write('# ' + str(int(iterations)) + ':\t' + str(maxL_new) + '\t' \
            + str(T) + '\t' + str(Prob(maxL, maxL_new, T)) + '\n')
        joints = joints_new
        maxL = maxL_new
    iterations = iterations + 1
f.close()


# ------ write the output ------

if html_output:
    OutputHTML(iterations)
    Draw(iterations)
else:
    OutputPlaintext()


print "analysis complete!"
