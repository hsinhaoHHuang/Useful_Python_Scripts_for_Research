#!/usr/bin/env python3

import unyt as u
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def Soliton_rho_c_from_r_c( m_22, r_c ):
    return ( 1.9e7*u.Msun/u.kpc**3 )*( m_22**-2 )*( (r_c/u.kpc)**-4 )

def Soliton_Density_profile( r, m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c( m_22, r_c )
    #return rho_c/( 1.0 + 9.1e-2*( r/r_c )**2 )**8
    return rho_c/( 1.0 + 9.050773267e-2*( r/r_c )**2 )**8

def Soliton_ShellDensity_profile( r, m_22, r_c ):
    return (4*np.pi*r**2)*Soliton_Density_profile( r, m_22, r_c )

def Soliton_EnclosedMass_profile( r, m_22, r_c ):
    M_enc = scipy.integrate.quad( lambda x: Soliton_ShellDensity_profile( x*r.units, m_22, r_c ), 0, r )[0]*( Soliton_ShellDensity_profile( r_c, m_22, r_c ).units*r_c.units )
    return M_enc

def Soliton_EnclosedMass_profile_1( r, m_22, r_c ):
    a           = ( 2.0**(1.0/8.0) - 1.0 )**(1.0/2.0)*( r/r_c )

    M_enc       = ( ( 4.2e9*u.Msun/( (10.0*m_22)**2 * (r_c/u.pc) ) )
                     *(1/(a**2 + 1)**7)
                     *(  3465*a**13 +
                        23100*a**11 +
                        65373*a**9 +
                       101376*a**7 +
                        92323*a**5 +
                        48580*a**3 -
                         3465*a +
                         3465*(a**2 + 1)**7*np.arctan(a)) )
    return M_enc

def Soliton_EnclosedMass_profile_2( r, m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c( m_22, r_c )
    gamma = r/r_c

    M_enc = 7.376*rho_c*(r_c**3)*( np.arctan(0.3017*gamma) + (1.0+0.091*gamma**2)**(-7)*
                                                                                        ( -3.017e-1*( gamma**(1)  ) +
                                                                                           3.849e-1*( gamma**(3)  ) +
                                                                                           6.656e-2*( gamma**(5)  ) +
                                                                                           6.651e-3*( gamma**(7)  ) +
                                                                                           3.903e-4*( gamma**(9)  ) +
                                                                                           1.255e-5*( gamma**(11) ) +
                                                                                           1.713e-7*( gamma**(13) )  ) )
    return M_enc

def Soliton_M_s_from_r_c_2( m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c( m_22, r_c )

    M_s   = 11.6*rho_c*(r_c**3)

    return M_s

def Soliton_rho_c_from_r_c_3( m_22, r_c ):
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar
    return (1.299/r_c)**4 / ( 4.0*np.pi*u.G*(ELBDM_Eta**2) )

def Soliton_M_s_from_r_c_3( m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c_3( m_22, r_c )

    M_s   = 25.9/(1.299**3)*rho_c*(r_c**3)
    return M_s

def Soliton_Ek_s_from_r_c_3( m_22, r_c ):
    M_s   = Soliton_M_s_from_r_c_3( m_22, r_c )
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    Ek_s  = M_s*(0.692/3.0*1.299**2)*(ELBDM_Eta*r_c)**-2
    return Ek_s

def Soliton_EkDensity_profile( r, m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c( m_22, r_c )
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    e_k   = 0.5 * (8.0*0.091)**2 * ELBDM_Eta**-2 * rho_c * r_c**-4 * r**2 /(1 +0.091*(r/r_c)**2 )**10
    return e_k

def Soliton_ShellEkDensity_profile( r, m_22, r_c ):
    return (4*np.pi*r**2)*Soliton_EkDensity_profile( r, m_22, r_c )

def Soliton_EnclosedEk_profile( r, m_22, r_c ):
    Ek_enc = scipy.integrate.quad( lambda x: Soliton_ShellEkDensity_profile( x*r.units, m_22, r_c ), 0, r )[0]*( Soliton_ShellEkDensity_profile( r_c, m_22, r_c ).units*r_c.units )
    return Ek_enc

def Soliton_sigma_from_r_c_3( m_22, r_c ):
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    sigma   = ( (2.0/9.0*0.692*1.299**2)*(ELBDM_Eta*r_c)**-2 )**0.5
    return sigma

def Soliton_M_c_from_CHrelation( m_22, M_h ):
    Mmin_0 = (4.4e7*u.Msun)*(m_22**-1.5)
    M_c = 0.25 * M_h**(1.0/3.0) * Mmin_0**(2.0/3.0)
    return M_c

def Soliton_r_c_from_CHrelation( m_22, M_h ):
    r_c = (1.6*u.kpc)/m_22*(M_h/(1.0e9*u.Msun))**(-1.0/3.0)
    return r_c

def Soliton_rho_c_from_CHrelation( m_22, M_h ):
    r_c   = Soliton_r_c_from_CHrelation( m_22, M_h )
    rho_c = Soliton_rho_c_from_r_c( r_c )
    return rho_c

def main() -> None:
    """main function
    This function is a template.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    m22   = 0.2
    Mh    = 1.88156694e+10*u.Msun
    rc    = Soliton_r_c_from_CHrelation( m22, Mh )
    rhoc  = Soliton_rho_c_from_r_c( m22, rc )
    rhoc3 = Soliton_rho_c_from_r_c_3( m22, rc )

    radius = np.linspace( 0.1*rc, 6.0*rc, 1000 )

    Menc   = u.unyt_array( [ Soliton_EnclosedMass_profile( r, m22, rc ) for r in radius.to('kpc') ] )
    Mc     = Soliton_EnclosedMass_profile( rc, m22, rc )
    Mc1    = Soliton_EnclosedMass_profile_1( rc, m22, rc )
    Mc2    = Soliton_EnclosedMass_profile_2( rc, m22, rc )
    Ms2    = Soliton_M_s_from_r_c_2( m22, rc )
    Ms3    = Soliton_M_s_from_r_c_3( m22, rc )
    Mch    = Soliton_M_c_from_CHrelation( m22, Mh )

    Ekenc  = u.unyt_array( [ Soliton_EnclosedEk_profile( r, m22, rc ) for r in radius.to('kpc') ] )
    Eks3   = Soliton_Ek_s_from_r_c_3( m22, rc )
    sigma3 = Soliton_sigma_from_r_c_3( m22, rc )

    # create figure
    fig = plt.figure(figsize=(10,7.5))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    with u.matplotlib_support:
        ax1.plot( radius.to('kpc'), Soliton_Density_profile( radius, m22, rc ).to('Msun/kpc**3'), linestyle='-',  linewidth=2.0, color='r', label='Soliton' )
        ax1.plot( [rc.to('kpc'), rc.to('kpc')], [np.min(Soliton_Density_profile( radius, m22, rc ).to('Msun/kpc**3')), np.max(Soliton_Density_profile( radius, m22, rc ).to('Msun/kpc**3')) ], linestyle='--',  linewidth=2.0, color='k', label='r_c = {: >4.3e}'.format( rc.to('kpc') ) )
        ax1.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [rhoc.to('Msun/kpc**3'),   rhoc.to('Msun/kpc**3')], linestyle='--',  linewidth=2.0, color='k', label=r'$\rho_c$ = {: >4.3e}'.format( rhoc.to('Msun/kpc**3') ) )
        ax1.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [rhoc3.to('Msun/kpc**3'), rhoc3.to('Msun/kpc**3')], linestyle='--',  linewidth=2.0, color='g', label=r'$\rho_c 3$ = {: >4.3e}'.format( rhoc3.to('Msun/kpc**3') ) )

        ax2.plot( radius.to('kpc'),                                              Menc.to('Msun'), linestyle='-',  linewidth=2.0, color='r', label='0' )
        ax2.plot( radius.to('kpc'), Soliton_EnclosedMass_profile_1( radius, m22, rc ).to('Msun'), linestyle='-',  linewidth=2.0, color='b', label='1' )
        ax2.plot( radius.to('kpc'), Soliton_EnclosedMass_profile_2( radius, m22, rc ).to('Msun'), linestyle='-',  linewidth=2.0, color='g', label='2' )
        ax2.plot( [rc.to('kpc'), rc.to('kpc')], [Soliton_EnclosedMass_profile( radius[0], m22, rc ).to('Msun'), Soliton_EnclosedMass_profile( radius[-1], m22, rc ).to('Msun') ], linestyle='--',  linewidth=2.0, color='k', label='r_c = {: >4.3e}'.format( rc.to('kpc') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Mc.to('Msun'),  Mc.to('Msun')],  linestyle='--',  linewidth=2.0, color='r', label=r'$M_c$ = {: >4.3e}'.format( Mc.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Mc1.to('Msun'), Mc1.to('Msun')], linestyle='--',  linewidth=2.0, color='b', label=r'$M_c, 1$ = {: >4.3e}'.format( Mc1.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Mc2.to('Msun'), Mc2.to('Msun')], linestyle='--',  linewidth=2.0, color='g', label=r'$M_c, 2$ = {: >4.3e}'.format( Mc2.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Ms2.to('Msun'), Ms2.to('Msun')], linestyle='-.',  linewidth=2.0, color='g', label=r'$M_s, 2$ = {: >4.3e}'.format( Ms2.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Ms3.to('Msun'), Ms3.to('Msun')], linestyle='-.',  linewidth=2.0, color='c', label=r'$M_s, 3$ = {: >4.3e}'.format( Ms3.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Mch.to('Msun'), Mch.to('Msun')], linestyle='--',  linewidth=2.0, color='m', label=r'$M_c, h$ = {: >4.3e}'.format( Mch.to('Msun') ) )

        ax3.plot( radius.to('kpc'), Ekenc.to('Msun*km**2/s**2'), linestyle='-',  linewidth=2.0, color='r', label='0' )
        ax3.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [Eks3.to('Msun*km**2/s**2'), Eks3.to('Msun*km**2/s**2')], linestyle='-.',  linewidth=2.0, color='c', label=r'$Ek_s, 3$ = {: >4.3e}'.format( Eks3.to('Msun*km**2/s**2') ) )
        ax3.plot( [rc.to('kpc'), rc.to('kpc')], [np.max(Ekenc.to('Msun*km**2/s**2')), np.min(Ekenc.to('Msun*km**2/s**2'))], linestyle='--',  linewidth=2.0, color='k', label='r_c = {: >4.3e}'.format( rc.to('kpc') ) )

        ax4.plot( radius.to('kpc'), np.sqrt( 2.0*Ekenc.to('Msun*km**2/s**2')/(3.0*Menc.to('Msun')) ), linestyle='-',  linewidth=2.0, color='r', label='0' )
        ax4.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [sigma3.to('km/s'), sigma3.to('km/s')], linestyle='-.',  linewidth=2.0, color='c', label=r'$\sigma, 3$ = {: >4.3e}'.format( sigma3.to('km/s') ) )
        ax4.plot( [rc.to('kpc'), rc.to('kpc')], [np.max(np.sqrt( 2.0*Ekenc.to('Msun*km**2/s**2')/(3.0*Menc.to('Msun')) )), np.min(np.sqrt( 2.0*Ekenc.to('Msun*km**2/s**2')/(3.0*Menc.to('Msun')) ))], linestyle='--',  linewidth=2.0, color='k', label='r_c = {: >4.3e}'.format( rc.to('kpc') ) )

    # x,y scale
    ax1.set_xscale( 'log' )
    ax1.set_yscale( 'log' )
    ax1.set_xlabel( r'$r \rm (kpc)$' )
    ax1.set_ylabel( r'$\rho$'+r' ($M_{\odot}/{\rm kpc}^3$)' )

    ax2.set_xscale( 'log' )
    ax2.set_yscale( 'log' )
    ax2.set_xlabel( r'$r \rm (kpc)$' )
    ax2.set_ylabel( r'$M$'+r' ($M_{\odot}$)' )

    ax3.set_xscale( 'log' )
    ax3.set_yscale( 'log' )
    ax3.set_xlabel( r'$r \rm (kpc)$' )
    ax3.set_ylabel( r'$E_k$'+r' ($M_{\odot}\rm km^2/s^2$)' )

    ax4.set_xscale( 'log' )
    ax4.set_yscale( 'log' )
    ax4.set_xlabel( r'$r \rm (kpc)$' )
    ax4.set_ylabel( r'$\sigma$'+r' ($\rm km/s$)' )

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()

    #ax.set_xlim( 0.0, 1.5*np.max(x) )

    # set the grid and ticks
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()

    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    ax1.tick_params( which='both',direction='in' )
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    ax2.tick_params( which='both',direction='in' )
    ax3.xaxis.set_ticks_position('both')
    ax3.yaxis.set_ticks_position('both')
    ax3.tick_params( which='both',direction='in' )
    ax4.xaxis.set_ticks_position('both')
    ax4.yaxis.set_ticks_position('both')
    ax4.tick_params( which='both',direction='in' )

    ## save to file
    plt.tight_layout( pad=0.1, h_pad=0.1, w_pad=0.1 )
    fig.savefig( 'fig_SolitonEnclosedMass.png' )

if __name__ == '__main__':
    main()
