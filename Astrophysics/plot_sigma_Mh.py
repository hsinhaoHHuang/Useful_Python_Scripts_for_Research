import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def power_law(x, a, b):
    return a*(x**b)

case    = np.array( [ ] )
m22     = np.array( [ ] )
Mh      = np.array( [ ] )
Rh      = np.array( [ ] )
vh      = np.array( [ ] )
sigma   = np.array( [ ] )
c_con   = np.array( [ ] )
c_fit   = np.array( [ ] )
r_gr    = np.array( [ ] )
r_gr_n  = np.array( [ ] )

marker_size = 60

fig = plt.figure(figsize=(8,6))

Sample_Mh_plot = np.logspace( 9, 12, num=1000 )

popt, pcov = curve_fit(power_law, Mh, sigma)

ax = fig.add_subplot(111)
ax.scatter( Mh, sigma, c=case, s=marker_size )
ax.plot( Sample_Mh_plot, 0.0100*(Sample_Mh_plot**(1.0/3.0) ), '--', label=r'$0.010 M_h^{1/3}$' )
ax.plot( Sample_Mh_plot, 0.0090*(Sample_Mh_plot**(1.0/3.0) ), '--', label=r'$0.009 M_h^{1/3}$' )
ax.plot( Sample_Mh_plot, 0.0080*(Sample_Mh_plot**(1.0/3.0) ), '--', label=r'$0.008 M_h^{1/3}$' )
ax.plot( Sample_Mh_plot, 0.0070*(Sample_Mh_plot**(1.0/3.0) ), '--', label=r'$0.007 M_h^{1/3}$' )
ax.plot( Sample_Mh_plot, popt[0]*(Sample_Mh_plot**(popt[1])), '--', label=r'$%.2f M_h^{%.2f}$'%(popt[0], popt[1]) )
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r'$M_h$'+" (Msun)")
ax.set_ylabel(r'$\sigma$'+" (km/s)")
ax.set_xlim( 1e9, 1e12 )
ax.set_ylim( 1e1, 1e2 )
ax.grid()
ax.legend()

fig.tight_layout( pad = 0.3, w_pad=-1.0 )
fig.savefig("fig_sigma_Mh.png")
plt.close()


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.scatter( c_fit, r_gr/Rh, c=case, s=marker_size )
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r'c')
ax.set_ylabel(r'$r_{\rm gr}/r_{\rm granule}$')
ax.grid()

fig.tight_layout( pad = 0.3, w_pad=-1.0 )
fig.savefig("fig_rgr_c.png")
plt.close()

