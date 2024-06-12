#!/usr/bin/env python3

import unyt as u
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def Soliton_rho_c_from_r_c_Schive2014( m_22, r_c ):
    return ( 1.9e7*u.Msun/u.kpc**3 )*( m_22**-2 )*( (r_c/u.kpc)**-4 )


def Soliton_rho_c_from_r_c_Chan2022( m_22, r_c ):
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    return (1.299/r_c)**4 / ( 4.0*np.pi*u.G*(ELBDM_Eta**2) )


def Soliton_Mass_Density_profile_Schive2014( r, m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c_Schive2014( m_22, r_c )

    return rho_c/( 1.0 + 9.1e-2*( r/r_c )**2 )**8


def Soliton_Mass_Enclosed_profile_NumInt( r, m_22, r_c ):

    M_enc = scipy.integrate.quad( lambda x: (4.0*np.pi*(x*r.units)**2)*Soliton_Mass_Density_profile_Schive2014( x*r.units, m_22, r_c ), 0, r )[0]*( Soliton_Mass_Density_profile_Schive2014( r_c, m_22, r_c ).units*r_c.units**3 )
    return M_enc


def Soliton_Mass_Enclosed_profile_Chiang2021( r, m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c_Schive2014( m_22, r_c )
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


def Soliton_Mass_Enclosed_profile_Chen2017( r, m_22, r_c ):
    a      = ( 2.0**(1.0/8.0) - 1.0 )**(1.0/2.0)*( r/r_c )

    M_enc  = ( ( 4.2e9*u.Msun/( (10.0*m_22)**2 * (r_c/u.pc) ) )
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


def Soliton_M_Total_from_r_c_Chiang2021( m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c_Schive2014( m_22, r_c )

    M_s   = 11.6*rho_c*(r_c**3)
    return M_s


def Soliton_M_Total_from_r_c_Chan2022( m_22, r_c ):
    rho_c = Soliton_rho_c_from_r_c_Chan2022( m_22, r_c )

    M_s   = 25.9/(1.299**3)*rho_c*(r_c**3)
    return M_s


def Soliton_Ek_Density_profile_AnaDif( r, m_22, r_c ):
    rho_c     = Soliton_rho_c_from_r_c_Schive2014( m_22, r_c )
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    e_k       = 0.5 * (8.0*0.091)**2 * ELBDM_Eta**-2 * rho_c * r_c**-4 * r**2 /(1 +0.091*(r/r_c)**2 )**10
    return e_k


def Soliton_Ek_Enclosed_profile_NumInt( r, m_22, r_c ):

    Ek_enc = scipy.integrate.quad( lambda x: (4*np.pi*(x*r.units)**2)*Soliton_Ek_Density_profile_AnaDif( x*r.units, m_22, r_c ), 0, r )[0]*( Soliton_Ek_Density_profile_AnaDif( r_c, m_22, r_c ).units*r_c.units**3 )
    return Ek_enc


def Soliton_Ek_Enclosed_profile_AnaInt( r, m_22, r_c ):
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar
    rho_c     = Soliton_rho_c_from_r_c_Schive2014( m_22, r_c )
    gamma     = r/r_c

    Ek_enc    = 2.0*np.pi*(8*0.091)**2 *(ELBDM_Eta)**-2 *rho_c*r_c*(   2.89556*( gamma**(17) ) +
                                                                     2.75767e2*( gamma**(15) ) +
                                                                     1.16088e4*( gamma**(13) ) +
                                                                     2.83024e5*( gamma**(11) ) +
                                                                     4.39244e6*( gamma**(9)  ) +
                                                                     4.48075e7*( gamma**(7)  ) +
                                                                     2.98080e8*( gamma**(5)  ) +
                                                                    -4.85617e8*( gamma**(3)  ) + ( 8.73479e-1*( gamma**(18) ) +
                                                                                                   8.63881e+1*( gamma**(16) ) +
                                                                                                   3.79728e+3*( gamma**(14) ) +
                                                                                                   9.73661e+4*( gamma**(12) ) +
                                                                                                   1.60494e+6*( gamma**(10) ) +
                                                                                                   1.76367e+7*( gamma**(8)  ) +
                                                                                                   1.29206e+8*( gamma**(6)  ) +
                                                                                                   6.08507e+8*( gamma**(4)  ) +
                                                                                                   1.67172e+9*( gamma**(2)  ) +
                                                                                                   2.04117e+9                   )*np.arctan(0.301662*gamma) - 6.15745e8*gamma )/(gamma**2 +10.989)**9
    return Ek_enc


def Soliton_Ek_Total_from_r_c_Chan2022( m_22, r_c ):
    M_s       = Soliton_M_Total_from_r_c_Chan2022( m_22, r_c )
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    Ek_s      = M_s*(0.692/3.0*1.299**2)*(ELBDM_Eta*r_c)**-2
    return Ek_s


def Soliton_sigma_from_r_c_Chan2022( m_22, r_c ):
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    sigma     = ( (2.0/9.0*0.692*1.299**2)*(ELBDM_Eta*r_c)**-2 )**0.5
    return sigma


def Soliton_Vintn_Shell_profile_AnaDif( r, m_22, r_c ):
    ELBDM_Eta = m_22*(1.0e-22*u.eV/u.c**2)/u.hbar

    Vintn = (8.0*0.091) * ELBDM_Eta**-1 * r_c**-2 * r /(1 +0.091*(r/r_c)**2 )
    return Vintn


def Soliton_M_c_from_CHrelation_Schive2014( m_22, M_h ):
    Mmin_0 = (4.4e7*u.Msun)*(m_22**-1.5)

    M_c    = 0.25*M_h**(1.0/3.0)*Mmin_0**(2.0/3.0)
    return M_c


def Soliton_r_c_from_CHrelation_Schive2014( m_22, M_h ):

    r_c = (1.6*u.kpc)/m_22*(M_h/(1.0e9*u.Msun))**(-1.0/3.0)
    return r_c


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

    m22    = 0.2
    Mh     = 1.88156694e+10*u.Msun
    rc     = Soliton_r_c_from_CHrelation_Schive2014( m22, Mh )

    radius = np.linspace( 0.1*rc, 4.0*rc, 1000 )

    # core density
    rhoc_Schive2014 = Soliton_rho_c_from_r_c_Schive2014( m22, rc )
    rhoc_Chan2022   = Soliton_rho_c_from_r_c_Chan2022( m22, rc )

    # core mass
    Menc_NumInt     = u.unyt_array( [ Soliton_Mass_Enclosed_profile_NumInt( r, m22, rc ) for r in radius.to('kpc') ] )
    Mc_NumInt       = Soliton_Mass_Enclosed_profile_NumInt(   rc, m22, rc )
    Mc_Chiang2021   = Soliton_Mass_Enclosed_profile_Chiang2021( rc, m22, rc )
    Mc_Chen2017     = Soliton_Mass_Enclosed_profile_Chen2017( rc, m22, rc )
    Ms_Chiang2021   = Soliton_M_Total_from_r_c_Chiang2021( m22, rc )
    Ms_Chan2022     = Soliton_M_Total_from_r_c_Chan2022( m22, rc )
    Mc_Schive2014   = Soliton_M_c_from_CHrelation_Schive2014( m22, Mh )

    # kinetic energy
    Ekenc_NumInt    = u.unyt_array( [ Soliton_Ek_Enclosed_profile_NumInt( r, m22, rc ) for r in radius.to('kpc') ] )
    Eks_Chan2022    = Soliton_Ek_Total_from_r_c_Chan2022( m22, rc )

    # velocity dispersion
    sigma_Chan2022  = Soliton_sigma_from_r_c_Chan2022( m22, rc )


    # create figure
    fig = plt.figure(figsize=(12,9))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    with u.matplotlib_support:
        ax1.plot(                                   radius.to('kpc'),                                                                                                 Soliton_Mass_Density_profile_Schive2014( radius, m22, rc ).to('Msun/kpc**3'),   linestyle='-',  linewidth=2.0, color='r', label='Schive 2014' )
        ax1.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                                  rhoc_Schive2014.to('Msun/kpc**3'),                                                    rhoc_Schive2014.to('Msun/kpc**3') ], linestyle='-.', linewidth=2.0, color='r', label=r'Schive 2014 $\rho_c$ = {: >4.3e}'.format( rhoc_Schive2014.to('Msun/kpc**3')  ) )
        ax1.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                                    rhoc_Chan2022.to('Msun/kpc**3'),                                                      rhoc_Chan2022.to('Msun/kpc**3') ], linestyle='-.', linewidth=2.0, color='b', label=r'Chan 2022   $\rho_c$ = {: >4.3e}'.format(   rhoc_Chan2022.to('Msun/kpc**3')    ) )
        ax1.plot( [           rc.to('kpc'),             rc.to('kpc')],  [np.min(Soliton_Mass_Density_profile_Schive2014( radius, m22, rc ).to('Msun/kpc**3')), np.max(Soliton_Mass_Density_profile_Schive2014( radius, m22, rc ).to('Msun/kpc**3'))], linestyle=':',  linewidth=2.0, color='k', label=r'$r_c$ = {: >4.3e}'.format( rc.to('kpc')                  ) )

        ax2.plot(                                   radius.to('kpc'),                                                   Menc_NumInt.to('Msun'),                                                                          linestyle='-',  linewidth=2.0, color='r', label='NumInt' )
        ax2.plot(                                   radius.to('kpc'),   Soliton_Mass_Enclosed_profile_Chiang2021( radius, m22, rc ).to('Msun'),                                                                          linestyle='--', linewidth=2.0, color='g', label='Chiang 2021' )
        ax2.plot(                                   radius.to('kpc'),     Soliton_Mass_Enclosed_profile_Chen2017( radius, m22, rc ).to('Msun'),                                                                          linestyle='--', linewidth=2.0, color='y', label='Chen 2017' )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                                 Mc_NumInt.to('Msun'),                                                   Mc_NumInt.to('Msun')], linestyle='-.', linewidth=2.0, color='r', label=r'NumInt      $M_c$ = {: >4.3e}'.format( Mc_NumInt.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                             Mc_Chiang2021.to('Msun'),                                               Mc_Chiang2021.to('Msun')], linestyle='-.', linewidth=2.0, color='g', label=r'Chiang 2021 $M_c$ = {: >4.3e}'.format( Mc_Chiang2021.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                               Mc_Chen2017.to('Msun'),                                                 Mc_Chen2017.to('Msun')], linestyle='-.', linewidth=2.0, color='y', label=r'Chen 2017   $M_c$ = {: >4.3e}'.format( Mc_Chen2017.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                             Mc_Schive2014.to('Msun'),                                               Mc_Schive2014.to('Msun')], linestyle='--', linewidth=2.0, color='m', label=r'CHrelation  $M_c$ = {: >4.3e}'.format( Mc_Schive2014.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                             Ms_Chiang2021.to('Msun'),                                               Ms_Chiang2021.to('Msun')], linestyle='-.', linewidth=2.0, color='g', label=r'Chiang 2021 $M_s$ = {: >4.3e}'.format( Ms_Chiang2021.to('Msun') ) )
        ax2.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [                                               Ms_Chan2022.to('Msun'),                                                 Ms_Chan2022.to('Msun')], linestyle='-.', linewidth=2.0, color='b', label=r'Chan 2022   $M_s$ = {: >4.3e}'.format( Ms_Chan2022.to('Msun') ) )
        ax2.plot( [           rc.to('kpc'),             rc.to('kpc')],  [ Soliton_Mass_Enclosed_profile_NumInt( radius[0], m22, rc).to('Msun'), Soliton_Mass_Enclosed_profile_NumInt( radius[-1], m22, rc ).to('Msun')], linestyle=':',  linewidth=2.0, color='k', label=r'$r_c$ = {: >4.3e}'.format(   rc.to('kpc') ) )

        ax3.plot(         radius.to('kpc'),                                                                                  Ekenc_NumInt.to('Msun*km**2/s**2'),   linestyle='-',  linewidth=2.0, color='r', label='NumInt' )
        ax3.plot(         radius.to('kpc'),                                         Soliton_Ek_Enclosed_profile_AnaInt( radius, m22, rc ).to('Msun*km**2/s**2'),   linestyle='--', linewidth=2.0, color='g', label='AnaInt' )
        ax3.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [        Eks_Chan2022.to('Msun*km**2/s**2'),         Eks_Chan2022.to('Msun*km**2/s**2') ], linestyle='-.', linewidth=2.0, color='b', label=r'Chan2022 $Ek_s$ = {: >4.3e}'.format( Eks_Chan2022.to('Msun*km**2/s**2') ) )
        ax3.plot( [           rc.to('kpc'),             rc.to('kpc')],  [ np.max(Ekenc_NumInt.to('Msun*km**2/s**2')), np.min(Ekenc_NumInt.to('Msun*km**2/s**2'))], linestyle=':',  linewidth=2.0, color='k', label=r'$r_c$ = {: >4.3e}'.format(                rc.to('kpc') ) )

        ax4.plot(         radius.to('kpc'),                               np.sqrt( 2.0*Ekenc_NumInt.to('Msun*km**2/s**2')/(3.0*Menc_NumInt.to('Msun')) ),                                                                                         linestyle='-',  linewidth=2.0, color='r', label='NumInt'           )
        ax4.plot(         radius.to('kpc'),                               np.sqrt( 2.0*Soliton_Ek_Enclosed_profile_AnaInt( radius, m22, rc ).to('Msun*km**2/s**2')/(3.0*Soliton_Mass_Enclosed_profile_Chiang2021( radius, m22, rc ).to('Msun'))), linestyle='--', linewidth=2.0, color='g', label='AnaInt'           )
        ax4.plot(         radius.to('kpc'),                               Soliton_Vintn_Shell_profile_AnaDif( radius, m22, rc ).to('km/s'),                                                                                                       linestyle='--', linewidth=2.0, color='c', label=r'$V_{\rm intn}$'  )
        ax4.plot( [np.min(radius.to('kpc')), np.max(radius.to('kpc'))], [ sigma_Chan2022.to('km/s'),                      sigma_Chan2022.to('km/s')],                                                                                             linestyle='-.', linewidth=2.0, color='b', label=r'Chan 2022 $\sigma$ = {: >4.3e}'.format( sigma_Chan2022.to('km/s') ) )
        ax4.plot( [rc.to('kpc'), rc.to('kpc')], [np.max(np.sqrt( 2.0*Ekenc_NumInt.to('Msun*km**2/s**2')/(3.0*Menc_NumInt.to('Msun')) )), np.min(np.sqrt( 2.0*Ekenc_NumInt.to('Msun*km**2/s**2')/(3.0*Menc_NumInt.to('Msun')) ))],                 linestyle=':',  linewidth=2.0, color='k', label=r'$r_c$ = {: >4.3e}'.format(        rc.to('kpc')  ) )

    # x,y scale
    #ax1.set_xscale( 'log' )
    #ax1.set_yscale( 'log' )
    ax1.set_xlabel( r'$r \rm (kpc)$' )
    ax1.set_ylabel( r'$\rho$'+r' ($M_{\odot}/{\rm kpc}^3$)' )

    # ax2.set_xscale( 'log' )
    # ax2.set_yscale( 'log' )
    ax2.set_xlabel( r'$r \rm (kpc)$' )
    ax2.set_ylabel( r'$M$'+r' ($M_{\odot}$)' )

    #ax3.set_xscale( 'log' )
    #ax3.set_yscale( 'log' )
    ax3.set_xlabel( r'$r \rm (kpc)$' )
    ax3.set_ylabel( r'$E_k$'+r' ($M_{\odot}\rm km^2/s^2$)' )

    #ax4.set_xscale( 'log' )
    #ax4.set_yscale( 'log' )
    ax4.set_xlabel( r'$r \rm (kpc)$' )
    ax4.set_ylabel( r'$\sigma$'+r' ($\rm km/s$)' )

    ax1.legend( loc='lower left'  )
    ax2.legend( loc='lower right' )
    ax3.legend( loc='lower right' )
    ax4.legend( loc='lower right' )

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
