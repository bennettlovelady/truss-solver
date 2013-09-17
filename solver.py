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

# uncomment this section to choose the output format/directory
'''
collating = True if raw_input('would you like to file these nicely? ') == 'y' else False

if collating:
    outpath = raw_input('input the output directory: ')
    if outpath[-1] != '/':
        outpath = outpath + '/'
    outputname = raw_input('what would you like to call these results? ')
else:
    outpath = inpath
    outputname = 'result'

html_output = True if raw_input('would you like html output? [y/n] ') == 'y' else False
'''


# ------ step 1 ------
print "reading input files.. "

with open(inpath+'joints.in','r') as fin:
    joints = np.genfromtxt(fin, comments="#", delimiter="\t")

with open(inpath+'members.in', 'r') as fin:
    members = np.genfromtxt(fin, comments="#", delimiter="\t")

with open(inpath+'supports.in', 'r') as fin:
    supports = np.genfromtxt(fin, comments="#")

with open(inpath+'applied.in','r') as fin:
    applied = np.genfromtxt(fin, comments="#", delimiter="\t")

material = {}
with open(inpath+'material.in', 'r') as fin:
    for line in fin:
        if (line[0] != '#'):
            spl = line.split()
            if len(spl)==2:
                (key, value) = spl 
                material[key] = float(value)


# ------ step 2 ------
print "constructing matrices.. "

zeroes = [[0]]*members.shape[0]
for n in range(5):
    members = np.concatenate((members, zeroes), 1)
# members table now has five extra columns (empty)

for k in range(members.shape[0]):
    i   = members[k,1]
    j   = members[k,2]
    dx  = joints[j-1,1] - joints[i-1,1]
    dy  = joints[j-1,2] - joints[i-1,2]
    L   = np.sqrt(dx**2 + dy**2)
    lij = dx/L
    mij = dy/L
    # put the values in the table
    members[k,3], members[k,4], members[k,5], members[k,6], members[k,7]\
        = dx, dy, L, lij, mij

# the members table is now fully populated with lij, mij


# ------ step 3: construct C ------
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
# C is complete!

'''
with open(inpath+'C.mat', 'w') as fout:
    np.savetxt(fout, C, fmt='%.3f', delimiter='\t')
'''


# ------ step 4: construct P ------
P = np.zeros(2*n)
num = applied.shape[0]-1
for i in range(1,num+1):
    P[2*applied[i,0]-2] = -applied[i,1]
    P[2*applied[i,0]-1] = -applied[i,2]


# ------ step 5: solve ------
print "solving.. "
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


# ------ step 6: find the maximum load ------
# use the material properties in materials.in to find the safe load
# go through every member and calculate its maximum safe load

maxL_t = np.zeros([members.shape[0],3])
for k in range(members.shape[0]):
    if Q[k] != 0:
        # yield: F/A <= s.f * s_y
        maxL_t[k,0] = material['sf'] * material['s_y'] * material['A'] / abs(Q[k])
        # buckling: F <= s.f * pi^2 * E * I / l^2
        if Q[k] < 0:
            maxL_t[k,1] = material['sf'] * np.pi**2 * material['E'] * material['I'] \
                            / (abs(Q[k]) * members[k,5]**2)
            # which is the smaller limit?
            maxL_t[k,2] = min(maxL_t[k,0], maxL_t[k,1])
        else:
            maxL_t[k,2] = maxL_t[k,0]

maxL = 1000000
maxL_m = []
for k in range(members.shape[0]):
    if maxL_t[k,2] != 0 and maxL_t[k,2] <= maxL:
        maxL = maxL_t[k,2]
for k in range(members.shape[0]):
    # to account for floating point errors,
    # check if the member's maxL is within 1% of maxL
    if np.allclose([maxL_t[k,2]], [maxL], rtol=1e-2):
        maxL_m.append(k)


# ------ step 7: present results ------
# 
# abandon all hope, ye who enter
#
# -----------------------------
print "writing output.. "

# some things used by both regular and html output sections:
# tension or compression?
def state(x):
    if (x<0):
        return 'comp.'
    elif (x>0):
        return 'tens.'
    else:
        return ''
    
# reaction force names
reactions = ['R'+str(int(supports[0]))+'x',\
             'R'+str(int(supports[1]))+'y',\
             'R'+str(int(supports[2]))+'y']
coord = ['x', 'y', 'y']

# maximum load of each member as strings
maxL_str = ['']*members.shape[0]
for k in range(members.shape[0]):
    maxL_str[k] = '---' if maxL_t[k,2] == 0 else str(round(maxL_t[k,2],2))

# maximum x-coord => bridge length
maxX = 0
for i in range(joints.shape[0]):
    if joints[i,1] > maxX:
        maxX = joints[i,1]

# total length of material used
sumL = 0
for k in range(members.shape[0]):
    sumL = sumL + members[k,5]

# regular, unformatted .out output
if not html_output:
    fout = open(outpath + outputname + '.out', 'w')
    if singular:
        fout.write("The system produced a singular matrix.\n\
	            These results were found by least squares approximation.\n")

    # output a table of all members' internal forces
    fout.write('mem. #\tstate\tforce/L\tmax load/N\n')
    for k in range(members.shape[0]):
        fout.write(str(k+1) + '\t' + state(Q[k]) + \
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
    fout.write('bridge length:\t' + str(int(maxX)) + 'mm\n')
    fout.write('total material:\t' + str(int(sumL)) + 'mm\n')

    fout.close()


# ------ step 8: draw the bridge ------
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
fig.savefig(outpath + outputname + '.png')


# ------ step 9: html page ------
# so far this is just a repeat of the normal output stage
def initTable(n):
    fout.write('<table>\n\t<colgroup>\n' + \
                n*'\t\t<col span="1" width="100">\n' + \
               '\t</colgroup>\n')

if html_output:
    fout = open(outpath + outputname + '.html', 'w')
    fout.write('<!DOCTYPE html>\n<html>\n<head>\n<title>Results</title>\n</head>\n')
    fout.write('\n<body>\n')
    fout.write('<img src="' + outputname + '.png">\n')
    if singular:
        fout.write('<p>The system produced a singular matrix.<br>\n' + \
	           'These results were found by least squares approximation.<br></p>\n\n')

    # table of internal forces
    fout.write('<p>Internal forces:</p>\n')
    initTable(4)
    fout.write('\t<tr><td>Member #</td><td>State</td>' + \
               '<td>Force/Load</td><td>Max Load/N</td></tr>\n')
    for k in range(members.shape[0]):
        fout.write('\t<tr><td>' + str(k+1) + \
                   '</td><td>' + state(Q[k]) + \
                   '</td><td>' + str(round(Q[k],4)) + \
                   '</td><td>' + maxL_str[k] + '</td></tr>\n')
    for i in range(3):
        fout.write('\t<tr><td>' + reactions[i] + \
                   '</td><td>' + '---' + \
                   '</td><td>' + str(round(Q[i-3],4)) + '</td><td></td></tr>\n')
    fout.write('</table>\n\n<br><br>')

    # weakest members
    fout.write('<p>Weakest members: </p>')
    initTable(3)
    fout.write('\t<tr><td>Member</td><td></td><td>Max Load / N</td></tr>\n')
    for i in range(len(maxL_m)):
        fout.write('\t<tr><td>Member ' + str(maxL_m[i]+1) + ':' + \
                   '</td><td>' + \
                   '</td><td>' + maxL_str[maxL_m[i]] + '</td></tr>\n')
    fout.write('</table><br><br>\n\n')

    # maximum load, bridge length and material used
    initTable(3)
    fout.write('\t<tr><td>Max load:</td><td></td><td>' + \
                str(round(maxL,4)) + ' N</td></tr>\n')
    fout.write('\t<tr><td></td><td></td><td>' + \
                str(round(maxL/9.81,4)) + ' kg</td></tr>\n')
    fout.write('\t<tr><td> </td></tr>\n')
    fout.write('\t<tr><td>Bridge length:</td><td></td><td>' + \
                str(int(maxX)) + 'mm</td></tr>\n')
    fout.write('\t<tr><td>Material used:</td><td></td><td>' + \
                str(int(sumL)) + 'mm</td></tr>\n')
    fout.write('</table><br><br>\n\n')

fout.close()

print "analysis complete!"
