import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from plot4d import plotter
from matplotlib.ticker import ScalarFormatter

case    = [ ]
m22     = [ ]
Mh      = [ ]
rho     = [ ]
vel     = [ ]
tau_gr  = [ ]
r_gr    = [ ]

marker_size = 60

fig = plt.figure(figsize=(9,4))

ax1 = fig.add_subplot(231)
ax1.scatter( m22, rho,  c=case, s=marker_size )
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel('$m_{22}$')
ax1.set_ylabel('$\\rho$'+" (Msun/kpc$^3$)")
ax1.xaxis.set_major_formatter(ScalarFormatter())
ax1.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)
ax1.set_xticks([0.2, 0.4, 0.8, 1.0])

ax2 = fig.add_subplot(232)
ax2.scatter( m22, vel,  c=case, s=marker_size )
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel('$m_{22}$')
ax2.set_ylabel('$v$'+" (km/s)")
ax2.set_xticks([0.2, 0.4, 0.8, 1.0])
ax2.xaxis.set_major_formatter(ScalarFormatter())
ax2.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)

ax3 = fig.add_subplot(233)
ax3.scatter( m22, Mh,  c=case, s=marker_size )
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlabel('$m_{22}$')
ax3.set_ylabel(r'$M_{\rm h}$'+" (Msun)")
ax3.set_xticks([0.2, 0.4, 0.8, 1.0])
ax3.xaxis.set_major_formatter(ScalarFormatter())
ax3.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)

ax4 = fig.add_subplot(234)
ax4.scatter( Mh, rho,  c=case, s=marker_size )
ax4.set_xscale("log")
ax4.set_yscale("log")
ax4.set_xlabel(r'$M_{\rm h}$'+" (Msun)")
ax4.set_ylabel('$\\rho$'+" (Msun/kpc$^3$)")

ax5 = fig.add_subplot(235)
ax5.scatter( Mh, rho,  c=case, s=marker_size )
ax5.set_xscale("log")
ax5.set_yscale("log")
ax5.set_xlabel(r'$M_{\rm h}$'+" (Msun)")
ax5.set_ylabel('$v$'+" (km/s)")

ax6 = fig.add_subplot(236)
ax6.scatter( rho, vel,  c=case, s=marker_size )
ax6.set_xscale("log")
ax6.set_yscale("log")
ax6.set_xlabel('$\\rho$'+" (Msun/kpc$^3$)")
ax6.set_ylabel('$v$'+" (km/s)")

fig.tight_layout( pad = 0.3, w_pad=-1.0 )
fig.savefig("fig_ParameterSpace_in.png")
plt.close()

#################################################################################################
fig = plt.figure(figsize=(9,4))

ax1 = fig.add_subplot(241)
ax1.scatter( m22, tau_gr,  c=case, s=marker_size )
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel('$m_{22}$')
ax1.set_ylabel(r'$\tau_{\rm gr}$'+" (Myr)")
ax1.set_xticks([0.2, 0.4, 0.8, 1.0])
ax1.xaxis.set_major_formatter(ScalarFormatter())
ax1.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)

ax2 = fig.add_subplot(242)
ax2.scatter( Mh, tau_gr,  c=case, s=marker_size )
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel(r'$M_{\rm h}$'+" (Msun)")
ax2.set_ylabel(r'$\tau_{\rm gr}$'+" (Myr)")
ax2.set_xticks([5e9, 1e10, 2e10])
ax2.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)

ax3 = fig.add_subplot(243)
ax3.scatter( rho, tau_gr,  c=case, s=marker_size )
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlabel('$\\rho$'+" (Msun/kpc$^3$)")
ax3.set_ylabel(r'$\tau_{\rm gr}$'+" (Myr)")

ax4 = fig.add_subplot(244)
ax4.scatter( vel, tau_gr,  c=case, s=marker_size )
ax4.set_xscale("log")
ax4.set_yscale("log")
ax4.set_xlabel('$v$'+" (km/s)")
ax4.set_ylabel(r'$\tau_{\rm gr}$'+" (Myr)")

ax5 = fig.add_subplot(245)
ax5.scatter( m22, r_gr,  c=case, s=marker_size )
ax5.set_xscale("log")
ax5.set_yscale("log")
ax5.set_xlabel('$m_{22}$')
ax5.set_ylabel(r'$r_{\rm gr}$'+" (kpc)")
ax5.set_xticks([0.2, 0.4, 0.8, 1.0])
ax5.xaxis.set_major_formatter(ScalarFormatter())
ax5.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)

ax6 = fig.add_subplot(246)
ax6.scatter( Mh, r_gr,  c=case, s=marker_size )
ax6.set_xscale("log")
ax6.set_yscale("log")
ax6.set_xlabel(r'$M_{\rm h}$'+" (Msun)")
ax6.set_ylabel(r'$r_{\rm gr}$'+" (kpc)")
ax6.set_xticks([5e9, 1e10, 2e10])
ax6.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)

ax7 = fig.add_subplot(247)
ax7.scatter( rho, r_gr,  c=case, s=marker_size )
ax7.set_xscale("log")
ax7.set_yscale("log")
ax7.set_xlabel('$\\rho$'+" (Msun/kpc$^3$)")
ax7.set_ylabel(r'$r_{\rm gr}$'+" (kpc)")

ax8 = fig.add_subplot(248)
ax8.scatter( vel, r_gr,  c=case, s=marker_size )
ax8.set_xscale("log")
ax8.set_yscale("log")
ax8.set_xlabel('$v$'+" (km/s)")
ax8.set_ylabel(r'$r_{\rm gr}$'+" (kpc)")

fig.tight_layout( pad = 0.3, w_pad=-1.0 )
fig.savefig("fig_ParameterSpace_out.png")
plt.close()
#################################################################################################

Log_Tau_gr = lambda x, y, z : np.log10(tau_gr[0]) + 3*(x-1-np.log10(m22[0])) + 6*(y+1-np.log10(vel[0])) -2*(z-np.log10(rho[0]))

gifname = plotter.plot4d( Log_Tau_gr, np.linspace(4,7,30), func_name="Log_Tau_gr", fps=5 )
