from pylab import *

### AUXILIARY FUNCTIONS
# 3PMMR curve for plotting
def fi_plot(x,a,b,c):
    f = ((-b-a*x)/c)**(-1)
    return f

# fix discontinuities in 3PMMRs for plotting
def fix_curve(yl):
    difs = np.roll(yl,-1)-yl
    difs = np.abs(difs)
    mask = np.where(difs>(.5*(np.max(yl)-np.min(yl))))
    yl2 = np.copy(yl)
    yl2[mask] = np.nan
    return(yl2)


# a lot of 3PMMRs (the ones tested for strength)
resl=['121','132','143','253','275','352','374','385','473','4117','3107','594','5149','3118','4139','7125','6115','495','583']
# our selected 3PMMRs
res_bl=['132','253','121','374','385','495','102','203','594']

l1,l2 = 1.2,1.7 # domain limits

# observed triplets
x0,y0 = np.genfromtxt('../datasets/triplets/known_lowmass_compact_triplets.dat',usecols=[0,1],delimiter=',',unpack=True)


### FIGURE
fac = 1.2
plt.figure(figsize=(10/fac,10/fac))

sl = 18
st = 15
se1 = 90
se2 = 50

xlabel(r'${\rm n_i/n_{(i+1)}}$',size=sl)
ylabel(r'${\rm n_{(i+1)}/n_{(i+2)}}$',size=sl)
xticks(size=st)
yticks(size=st)

# observed triplets
plt.scatter(x0,y0,s=se1,color='yellow',zorder=112,label='Observados')
plt.scatter(x0,y0,s=se2,color='black',zorder=112,label='Observados')

# plot 2PMMR
for res2p in [5/4.,4/3.,3/2.]:
    plt.vlines(res2p,l1,l2,linestyles='solid',colors='black',lw=1)
    plt.hlines(res2p,l1,l2,linestyles='solid',colors='black',lw=1)

# 3PMMRs - strength width and plotting
dom=np.linspace(l1,l2+.2,50)
for res in resl:
    # samples of the strength F of the 3PMMR (a,b,c) at point (x,y) of the period ratio plane
    x,y,a,b,c,F=np.genfromtxt('mmrs_sampled/'+res+'.dat',unpack=True,skip_header=1) 
    a,b,c=int(a[0]),int(b[0]),int(c[0])

    inplot = (x>=l1/2.)&(x<=l2*2)&(y>=l1/2.)&(y<=l2*2)

    y = y[inplot]
    x = x[inplot]
    F = F[inplot]
    
    # width proportional to F by an arbitrary function that allows for visual comparison
    size=2.5*(F*10**9)*2.  #
    size*=0.2
    
    # plot widths
    if res in res_bl:
        color = '#86b5d5'
    else:
        color = 'lightcoral'
    plt.fill_between(x, y-size*.5, y+size*.5,alpha=1,color=color,edgecolor=None)
    plt.fill_betweenx(y,x-size*.5, x+size*.5,alpha=1,color=color,edgecolor=None)

    # grey areas
    if res=='4139':
        ylimmin = y-0.01-0.04*(y-1.33)
        #        plot(x,ylimmin,color='black',ls='dashed',zorder=1000)
        fill_between(x,ylimmin,color='gray',alpha=0.3,zorder=300)

    if res=='594':
        ylimmax = y+0.02+0.04*(y-1.5)
        #        plot(x,ylimmax,color='black',ls='dashed',zorder=1000)
        fill_between(x,ylimmax,1.7,color='gray',alpha=0.3,zorder=300)
    
    # MMR labels
    a=int(res[0])
    b=-1*int(res[1])
    c=int(res[-1])
    if len(res)==4:
        b=-1*int(res[1:3])        
    
    fi = fi_plot(dom,a,b,c) # y values of (a,b,c)
    xtext = dom[(fi<=l2)&(fi>=l1)&(dom<=l2)&(dom>=l1)][-1]+0.005
    ytext = fi_plot(xtext,a,b,c)
    
    text = '('+str(a)+','+str(b)+','+str(c)+')'
    
    if ytext>=1.66:
        xtext = xtext - 0.015
        ytext = l2 + .005
        text = '('+str(a)+','+str(b)+','+str(c)+')'
    if res=='473' or res=='121':
        xtext+=0.01
    if res=='4117':
        xtext = dom[(fi<=l2)&(fi>=l1)&(dom<=l2)&(dom>=l1)][-1]+0.004
        ytext = fi_plot(xtext,a,b,c)
        text = '('+str(a)+','+str(b)+','+str(c)+')'
    
    if res=='7125' or res=='6115' or res=='4117' or res=='4139' or res=='5149' or res=='3107' or res=='3118' or res=='583' or res=='352' or res=='473' or res=='275' or res=='143' or res=='':
        text = ' '
        
    rot = 60
    if res == '583' or res == '352':
        rot = 80
    if res == '473' or res == '594':
        rot = 76
    if res == '121' :
        rot = 72
    if res == '374' or res == '495' :
        rot = 65
    if res == '132' or res == '275' :
        rot = 45
    if res == '275' :
        rot = 37
    if res == '143' :
        rot = 30
        
    plt.text(xtext,ytext,text,size=14,rotation=rot)

plt.xlim(l1,l2)
plt.ylim(l1,l2)

savefig('res_strengths.png',bbox_inches='tight',dpi=100)

show()
