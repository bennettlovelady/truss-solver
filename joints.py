#! /usr/bin/env python

print "this script converts integer joint coordinates to mm"
print "file format is [joint #] TAB [xpos] TAB [ypos]"
print "but just use spaces here"
print "input -1 -1 as the x,y coords to signify end of input"
path = raw_input("what is the working directory to place joints.in: ")

scale = int(raw_input("enter the unit length in mm: "))

if path[-1]!='/':
    path = path + '/'
f = open(path+'joints.in', 'w')

n = 1
p = [0,0]
while p[0]!='-1': 
    p = raw_input(str(n) + ': ').split()
    if p[0]!='-1':
        f.write(str(n) + '\t' \
                       + str(int(scale*float(p[0]))) + '\t' \
                       + str(int(scale*float(p[1]))) + '\n')
        n = n+1

f.close()
