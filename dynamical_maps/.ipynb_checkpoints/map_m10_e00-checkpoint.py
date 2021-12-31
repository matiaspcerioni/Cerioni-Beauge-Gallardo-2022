import pylab
from pylab import *
import numpy as np

# READ DATA
mapa = genfromtxt("simulations/chaos_tot_m10_e00.dat")
data = genfromtxt('../datasets/triplets/known_lowmass_compact_triplets.dat',usecols=[0,1],delimiter=',',unpack=True)

# PARAMETERS
g   = 0.1720209895**2
m0  = 1.0
m   = 0.0001 # solar masses
mu1 = g*(m0+m)
mu2 = g*(m0+m+m)
mu3 = g*(m0+m+m+m)

nx = 1000
ny = nx

# DYNAMICAL MAP
a1 = (mapa[:,3][mapa[:,1]==1]).reshape(ny,nx)
a2 = (mapa[:,3][mapa[:,1]==2]).reshape(ny,nx)
a3 = (mapa[:,3][mapa[:,1]==3]).reshape(ny,nx)
delta_a1 = (mapa[:,9][mapa[:,1]==1]).reshape(ny,nx)
delta_a2 = (mapa[:,9][mapa[:,1]==2]).reshape(ny,nx)
delta_a3 = (mapa[:,9][mapa[:,1]==3]).reshape(ny,nx)
delta_e1 = (mapa[:,10][mapa[:,1]==1]).reshape(ny,nx)
delta_e2 = (mapa[:,10][mapa[:,1]==2]).reshape(ny,nx)
delta_e3 = (mapa[:,10][mapa[:,1]==3]).reshape(ny,nx)
delta_e = np.ndarray([ny,nx])
delta_a = np.ndarray([ny,nx])
for i in range(nx):
    for j in range(ny):
        delta_e[i,j] = max(delta_e1[i,j],delta_e2[i,j],delta_e3[i,j])
        delta_a[i,j] = max(delta_a1[i,j],delta_a2[i,j],delta_a3[i,j])
p1p2 = sqrt(((1.0+2.0*m)/(1.0+1.0*m))*((a1/a2)**3))
p2p3 = sqrt(((1.0+3.0*m)/(1.0+2.0*m))*((a2/a3)**3))
n1n2 = 1.0/p1p2
n2n3 = 1.0/p2p3
p12 = 1.0/(1.0/p1p2 - 1.0)
p13 = 1.0/(1.0/p1p2/p2p3 - 1.0)
p23 = 1.0/(1.0/p2p3 - 1.0)
n1 = sqrt(mu1/a1/a1/a1)
n2 = sqrt(mu2/a2/a2/a2)
n3 = sqrt(mu2/a3/a3/a3)
del_l1 = g*m*m*(1.43*p12+0.13)/(n1-n2)/a2 + g*m*m*(1.43*p13+0.13)/(n1-n3)/a3
del_l2 = g*m*m*(1.43*p12+0.13)/(n1-n2)/a2 + g*m*m*(1.43*p23+0.13)/(n2-n3)/a3
del_l3 =-g*m*m*(1.43*p23+0.13)/(n2-n3)/a3 - g*m*m*(1.43*p13+0.13)/(n1-n3)/a3
del_a1 = 2.0*sqrt(a1/mu1)*abs(del_l1)/m
del_a2 = 2.0*sqrt(a2/mu2)*abs(del_l2)/m
del_a3 = 2.0*sqrt(a3/mu3)*abs(del_l3)/m
del_a = np.ndarray([ny,nx])
for i in range(nx):
    for j in range(ny):
        del_a[i,j] = max(del_a1[i,j],del_a2[i,j],del_a3[i,j])

# FIGURE
fac = 1.2
fig = figure(figsize=(10/fac,10/fac))
#
xymin = 1.2
xymax = 1.7
#
vmin = -8.0
vmax = -3.0
#
ecol = 'yellow'
sl = 20
st = 16
se1 = 90
se2 = 50
#
ax = subplot(1,1,1)
xlim(xymin,xymax)
ylim(xymin,xymax)
xlabel('n$_1$/n$_2$',size=sl)
ylabel('n$_2$/n$_3$',size=sl)
xticks(size=st)
yticks(size=st)
#
pcolormesh(n1n2,n2n3,log(abs(delta_a)),vmin=vmin,vmax=vmax,cmap='twilight_shifted',shading='gouraud')
#
scatter(data[0],data[1],s=se1,color=ecol)
scatter(data[0],data[1],s=se2,color='black')
#
tight_layout(pad=0.5,h_pad=1)
savefig("map_m10_e00.png",bbox_inches='tight',dpi=150)
