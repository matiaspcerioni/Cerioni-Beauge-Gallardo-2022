import pylab
from pylab import *
import numpy as np

# READ DATA
mapa05 = genfromtxt("simulations/chaos_tot_m10_e05.dat")
data = genfromtxt('../datasets/triplets/known_lowmass_compact_triplets.dat',usecols=[0,1],delimiter=',',unpack=True)

# PARAMETERS
g   = 0.1720209895**2
m0  = 1.0
m05   = 10.0*3.0404e-6 # solar masses
mu051 = g*(m0+m)
mu052 = g*(m0+m+m)
mu053 = g*(m0+m+m+m)

nx = 1000
ny = nx

# DYNAMICAL MAP
a051 = (mapa05[:,3][mapa05[:,1]==1]).reshape(ny,nx)
a052 = (mapa05[:,3][mapa05[:,1]==2]).reshape(ny,nx)
a053 = (mapa05[:,3][mapa05[:,1]==3]).reshape(ny,nx)
delta_a105 = (mapa05[:,9][mapa05[:,1]==1]).reshape(ny,nx)
delta_a205 = (mapa05[:,9][mapa05[:,1]==2]).reshape(ny,nx)
delta_a305 = (mapa05[:,9][mapa05[:,1]==3]).reshape(ny,nx)
delta_e105 = (mapa05[:,10][mapa05[:,1]==1]).reshape(ny,nx)
delta_e205 = (mapa05[:,10][mapa05[:,1]==2]).reshape(ny,nx)
delta_e305 = (mapa05[:,10][mapa05[:,1]==3]).reshape(ny,nx)
delta_e05 = np.ndarray([ny,nx])
delta_a05 = np.ndarray([ny,nx])
for i in range(nx):
    for j in range(ny):
        delta_e05[i,j] = max(delta_e105[i,j],delta_e205[i,j],delta_e305[i,j])
        delta_a05[i,j] = max(delta_a105[i,j],delta_a205[i,j],delta_a305[i,j])
p1p205 = sqrt(((1.0+2.0*m05)/(1.0+1.0*m05))*((a051/a052)**3))
p2p305 = sqrt(((1.0+3.0*m05)/(1.0+2.0*m05))*((a052/a053)**3))
n1n205 = 1.0/p1p205
n2n305 = 1.0/p2p305
p0512 = 1.0/(1.0/p1p205 - 1.0)
p0513 = 1.0/(1.0/p1p205/p2p305 - 1.0)
p0523 = 1.0/(1.0/p2p305 - 1.0)
n051 = sqrt(mu051/a051/a051/a051)
n052 = sqrt(mu052/a052/a052/a052)
n053 = sqrt(mu052/a053/a053/a053)
del05_l1 = g*m05*m05*(1.43*p0512+0.13)/(n051-n052)/a052 + g*m05*m05*(1.43*p0513+0.13)/(n051-n053)/a053
del05_l2 = g*m05*m05*(1.43*p0512+0.13)/(n051-n052)/a052 + g*m05*m05*(1.43*p0523+0.13)/(n052-n053)/a053
del05_l3 =-g*m05*m05*(1.43*p0523+0.13)/(n052-n053)/a053 - g*m05*m05*(1.43*p0513+0.13)/(n051-n053)/a053
del05_a1 = 2.0*sqrt(a051/mu051)*abs(del05_l1)/m05
del05_a2 = 2.0*sqrt(a052/mu052)*abs(del05_l2)/m05
del05_a3 = 2.0*sqrt(a053/mu053)*abs(del05_l3)/m05
del05_a = np.ndarray([ny,nx])
for i in range(nx):
    for j in range(ny):
        del05_a[i,j] = max(del05_a1[i,j],del05_a2[i,j],del05_a3[i,j])

# FIGURE
fac = 1.2
fig = figure(figsize=(10/fac,10/fac))
#
ecol = 'yellow'
sl = 20
st = 16
#
ax = subplot(1,1,1)
xymin = 1.2
xymax = 1.7
xlim(xymin,xymax)
ylim(xymin,xymax)
xlabel('n$_1$/n$_2$',size=sl)
ylabel('n$_2$/n$_3$',size=sl)
xticks(size=st)
yticks(size=st)
#
vmin = -8.0
vmax = -3.5
pcolormesh(n1n205,n2n305,log(abs(delta_a05)),vmin=vmin,vmax=vmax,cmap='twilight_shifted',shading='gouraud')#
#
se1 = 90
se2 = 50
scatter(data[0],data[1],s=se1,color=ecol)
scatter(data[0],data[1],s=se2,color='black')
#
tight_layout(pad=0.5,h_pad=1)
#
savefig("map_m10_e05.png",bbox_inches='tight',dpi=150)
