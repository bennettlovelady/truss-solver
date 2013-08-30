#! /usr/bin/env python

'''
solves a truss defined by the files joints.in, members.in and supports.in
returns a file forces.out with the resultant forces

step 1: read input files joints.in, members.in, supports.in, applied.in.
        begin constructing the relevant tables
step 2: expand the members table to include dx, dy, L, lij, mij
step 3: construct a 2nx2n matrix of coefficients C, where n = # joints
step 4: construct the "applied forces" vector P
step 5: solve the matrix equation C.Q=P, where Q is a list of unknowns
'''

import numpy as np


# ------ step 1 ------
fin = open('joints.in', 'r')
joints = np.genfromtxt(fin, comments="#", delimiter="\t")
fin.close()

fin = open('members.in', 'r')
members = np.genfromtxt(fin, comments="#", delimiter="\t")
fin.close()

fin = open('supports.in', 'r')
supports = np.genfromtxt(fin, comments="#")
fin.close()

fin = open('applied.in','r')
applied = np.genfromtxt(fin, comments="#", delimiter="\t")
fin.close()


# ------ step 2 ------
zeroes = [[0]]*members.shape[0]
for n in range(5):
    members = np.concatenate((members, zeroes), 1)

# members table now has five extra columns (empty)

for n in range(members.shape[0]):
    i   = members[n,1]
    j   = members[n,2]
    # get the values
    dx  = joints[j-1,1] - joints[i-1,1]
    dy  = joints[j-1,2] - joints[i-1,2]
    L   = np.sqrt(dx**2 + dy**2)
    lij = dx/L
    mij = dy/L
    # put them in the table
    members[n,3], members[n,4], members[n,5], members[n,6], members[n,7]\
        = dx, dy, L, lij, mij

'''
fout = open('out.out','w')
fout.write(np.array_str(members, max_line_width=150))
fout.close()
'''
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
C[2*(supports[0])-2, 2*n-3] = 1
C[2*(supports[1])-1, 2*n-2] = 1
C[2*(supports[2])-1, 2*n-1] = 1
'''
fout = open('c.out', 'w')
fout.write(np.array_str(C, max_line_width=500))
fout.close()
'''
# C is complete!


# ------ step 4: construct P ------
P = np.zeros(2*n)
num = applied.shape[0]-1
for i in range(1,num+1):
    P[2*applied[i,0]-2] = -applied[i,1]
    P[2*applied[i,0]-1] = -applied[i,2]


# ------ step 5: solve ------
Q = np.linalg.solve(C,P)

names = ['R'+str(int(supports[0]))+'x',\
         'R'+str(int(supports[1]))+'y',\
         'R'+str(int(supports[2]))+'y']
coord = ['x', 'y', 'y']

def state(x):
    if (x<0):
        return 'comp.'
    elif (x>0):
        return 'tens.'
    else:
        return ''

fout = open('result.out', 'w')

for i in range(members.shape[0]):
    fout.write(str(i+1) + '\t' + state(Q[i]) + '\t' + str(Q[i]) + '\n')
for i in range(3):
    x = members.shape[0] + i
    fout.write('R' + str(int(supports[i])) + coord[i] + '\t' \
                   + state(Q[x]) + '\t' \
                   + str(Q[x]) + '\n')

fout.close()

