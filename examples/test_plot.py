


import matplotlib as mpl
import numpy as np
import flare.plt as fplt



# --- line, fill_between

x = np.array([7,12])
y = np.array([28,31])

ax.fill_between(x,y-0.2,y+0.2, color='k', alpha=0.2)
ax.plot(x,y, c='k', lw=2)

# ---

for m in np.linspace(0.0,0.5,5):

    ax.plot(x,(x-7)*m+27.5, c='k', lw=1)



# --- scatter1

N = 50
x = 8 + np.random.randn(N)/5
y = 30 + np.random.randn(N)/5
z = x**2 + y**2

ax.scatter(x,y,c=z,cmap=cm.plasma,s=10)


ax.set_ylim([27, 32])
ax.set_xlim([7, 12])

ax.set_ylabel(r'$\log_{10}(L_{\nu}/{\rm erg\ s^{-1}\ Hz^{-1})}$')
ax.set_xlabel(r'$\log_{10}(M_{\star}/{\rm M_{\odot}})$')

ax.grid(True)




fig.savefig('test_plot.pdf')
