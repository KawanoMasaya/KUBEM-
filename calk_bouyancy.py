import linecache
import numpy as np
import math

NBODY=input('input number of body:')

for IBODY in range(int(NBODY)):

	print('input filename of the model (.dat)')
	for i in range(10):
		infil0=input()
		if '#' in infil0:
			continue
		break

	infil1 = infil0+str('.dat')
	open(infil1)

	print('input file name for output (.out)')
	for i in range(10):
		filen0 = input()
		if '#' in filen0:
			continue
		break

	filen1 = filen0+str('.out')
	open(filen1 , 'w')

	open('output.out','w')

	XLOC,YLOC=input('input XLOC , YLOC').split()


print('\nremove Irregular Frequency (yes=1,no=0)')
for i in range(10):
	IR=input()
	if '#' in IR:
		continue
	break
if int(IR) !=0 and int(IR) !=1:
	print('parameter error for IR')
	exit

print('output data file (yes=1,no=0,single line)')

for i in range(30):
	expt=input()
	if '#' in expt:
		continue
	break
expt=expt.split()

for i in range(12):
	expt[i]=int(expt[i])

for i in range(12):
	if int(expt[i]) !=0 and int(expt[i]) !=1:
		print('parameter error for expt')
		exit



NINP=open(str(infil1))
for a in NINP:
	s=NINP.readline()
	if '#' in s:
		continue
	break
title=s
out=open(str(filen1),'a')
out.write(title+'\n')

out.write('GENERAL PARAMETERS'+'\n')


for a in NINP:
	s=NINP.readline()
	if '#' in s:
		continue
	break
nsym,node0,ne0,nodei0,nei0,wp,rho,g = s.split()
nsym = int(nsym)
node0 = int(node0)
ne0 = int(ne0)
nodei0 = int(nodei0)
nei0 = int(nei0)
wp = float(wp)
rho = float(rho)
g = float(g)

for a in NINP:
	s=NINP.readline()
	if '#' in s:
		continue
	break
Xoffset,Yoffset,Zoffset,kxx,kyy,kzz = s.split()

mode = 6*int(NBODY)

if int(nsym) == 1:
	out.write(' symmetry         =   ASYMMETRY'+'\n')
elif int(nsym) == 2:
	out.write(' symmetry         =          YZ'+'\n')
elif int(nsym) == 4:
	out.write(' symmetry         =       XZ,YZ'+'\n')
elif int(nsym) ==3:
	print('nsym = 3(XZ) is in preparation')
	exit

out.write(' Irregular Freq   =           '+str(IR)+'\n')
out.write(' node number(0)   =           '+str(node0)+'\n')
out.write(' element number(0)=           '+str(ne0)+'\n')
out.write(' node number(i)   =           '+str(nodei0)+'\n')
out.write(' element number(i)=           '+str(nei0)+'\n')
if float(wp) <= 0:
	out.write(' water depth(m)   =          Infinity'+'\n')
else:
	out.write(' water depth(m)   =           '+str(wp)+'\n')
out.write(' rho              =           '+str(rho)+'\n')
out.write(' g                =           '+str(g)+'\n')
out.write('position of center of gravity'+'\n')
out.write(' Gx               =           '+str(Xoffset)+'\n')
out.write(' Gy               =           '+str(Yoffset)+'\n')
out.write(' Gz               =           '+str(Zoffset)+'\n')
out.write('gyrational radius'+'\n')
out.write(' kxx              =           '+str(kxx)+'\n')
out.write(' kyy              =           '+str(kyy)+'\n')
out.write(' kzz              =           '+str(kzz)+'\n')




#BMAIN

if nsym != 3:
	ne =ne0*nsym
	nei = nei0*nsym
else:
	ne = ne0*(nsym-1)
	nei = nei0*(nsym-1)

if expt[1] ==1:
	open('01_mesh.out', 'w')
	meshf=open('01_mesh.out','a')
	meshf.write('BODY SURFACE'+'\n')

#ELMINP
l = 0
lg = 0
xi=[0,0,0]
X = np.zeros((3,node0))
Rt = [0]*node0

for i in range(100000):
	s=NINP.readline()
	if '#' in s:
		continue
	n,ng,xi[0],xi[1],xi[2],rtz=s.split()

	n = int(n)
	ng = int(ng)
	xi[0] = float(xi[0])
	xi[1] = float(xi[1])
	xi[2] = float(xi[2])
	rtz = int(rtz)

	X[0,n-1] = xi[0]
	X[1,n-1] = xi[1]
	X[2,n-1] = xi[2]
	Rt[n-1] = rtz


	if lg != 0:
		lg = abs(lg)*np.sign(n-1)
		li = (abs(n-1+lg)-1)/abs(lg)
		xi[0] = (X[0,n-1]-X[0,l-1])/float(li)
		xi[1] = (X[1,n-1]-X[1,l-1])/float(li)
		xi[2] = (X[2,n-1]-X[2,l-1])/float(li)

		for k in range(100000):
			l = l+lg
			if (n-1)*lg <= 0:
				break
			if l <= 0 or l > node0:
				print('error at rdnode')
				exit
			X[0,l-1] = X[0,l-lg-1] + xi[0]
			X[1,l-1] = X[1,l-lg-1] + xi[1]
			X[2,l-1] = X[2,l-lg-1] + xi[2]
			Rt[l-1] = Rt[l-lg-1]

	if n == node0:
		break

	l = n
	lg = ng

for i in range(node0):
	if Rt[i-1] == 1:
		r = X[0,i-1]
		if X[1,i-1] != 0:
			th = [1,i-1]*math.pi/180
			X[0,i-1] = r*math.cos(th)
			X[1,i-1] = r*math.sin(th)
		else:
			X[0,i-1] = 0
			X[1,i-1] = r

for i in range(node0):
	X[0,i-1] = X[0,i-1] + float(Xoffset)
	X[1,i-1] = X[1,i-1] + float(Yoffset)

if expt[1] == 1:
	meshf.write(' coodinate of nodes')

	for k in range(node0):
		print(k)
