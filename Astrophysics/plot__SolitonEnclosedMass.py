
import numpy as np
import scipy.integrate

m_22   = 0.2
M_h    = 1.88156694e+10
Mmin_0 = 4.91991577e+08

def Soliton_rho_c_from_r_c( r_c ):   # in Msun/kpc^3
    return 1.9e7/( m_22**2 * r_c**4 )

def Soliton_M_c_from_r_c( r_c, r ):     # in Msun
    constfactor = ( 2.0**(1.0/8.0) - 1.0 )**(1.0/2.0)*( r/r_c )
    M_c         = ( ( 4.2e9/( (10.0*m_22)**2 * (r_c*1.0e3) ) )
                     *(1/(constfactor**2 + 1)**7)
                     *(3465*constfactor**13 + 23100*constfactor**11 + 65373*constfactor**9 + 101376*constfactor**7 + 92323*constfactor**5 + 48580*constfactor**3 - 3465*constfactor + 3465*(constfactor**2 + 1)**7*np.arctan(constfactor)) )
    return M_c

def Soliton_r_c_M_c_from_r_c( r_c ): # in Msun*kpc
    return r_c*Soliton_M_c_from_r_c( r_c, r_c )

def Soliton_M_c_from_r_c_Chiang( r_c, r ):
    rho_c = Soliton_rho_c_from_r_c( r_c )
    gamma = r/r_c

    M_c = 7.376*rho_c*(r_c**3)*( np.arctan(0.3017*gamma) + (1.0+0.091*gamma**2)**(-7)*
                                                                                      ( -3.017e-1*( gamma**(1)  ) +
                                                                                         3.849e-1*( gamma**(3)  ) +
                                                                                         6.656e-2*( gamma**(5)  ) +
                                                                                         6.651e-3*( gamma**(7)  ) +
                                                                                         3.903e-4*( gamma**(9)  ) +
                                                                                         1.255e-5*( gamma**(11) ) +
                                                                                         1.713e-7*( gamma**(13) )  ) )
    return M_c

def Soliton_r_c_M_c_from_r_c_Chiang( r_c ): # in Msun*kpc
    return r_c*Soliton_M_c_from_r_c_Chiang( r_c, r_c )

def Soliton_M_s_from_r_c_Chiang( r_c ):
    rho_c = Soliton_rho_c_from_r_c( r_c )

    M_s   = 11.6*rho_c*(r_c**3)

    return M_s

def Soliton_r_c_M_s_from_r_c_Chiang( r_c ): # in Msun*kpc
    return r_c*Soliton_M_s_from_r_c_Chiang( r_c )

def Soliton_M_s_from_r_c_Chan( r_c ):
    rho_c = Soliton_rho_c_from_r_c( r_c )

    M_s   = 25.9/(1.299**3)*rho_c*(r_c**3)

    return M_s

def Soliton_r_c_M_s_from_r_c_Chan( r_c ): # in Msun*kpc
    return r_c*Soliton_M_s_from_r_c_Chan( r_c )

def Soliton_r_c_M_s_from_r_c_Chan_2( r_c ):
    Newton_G  = 4.51696776e-39
    ELBDM_Eta = 3.21902202e+14
    return 25.9*1.299/( 4.0*np.pi*Newton_G*(ELBDM_Eta**2) )

def Soliton_rho_c_from_r_c_Chan( r_c ):
    Newton_G  = 4.51696776e-39
    ELBDM_Eta = 3.21902202e+14
    return (1.299/r_c)**4 / ( 4.0*np.pi*Newton_G*(ELBDM_Eta**2) )

def Soliton_analytical_dens( r, r_c ):  # in Msun/kpc^3
    rho_c = Soliton_rho_c_from_r_c( r_c )
    #return rho_c/( 1.0 + 9.1e-2*( r/r_c )**2 )**8
    #return rho_c/( 1.0 + 9.051e-2*( r/r_c )**2 )**8
    return rho_c/( 1.0 + 9.050773267e-2*( r/r_c )**2 )**8

def Soliton_analytical_shell_density( r, r_c ):
    return (4*np.pi*r**2)*Soliton_analytical_dens( r, r_c )

def Soliton_M_c_from_r_c_NumInt( r, r_c ):
    M_c = scipy.integrate.quad( lambda x: Soliton_analytical_shell_density( x, r_c ), 0, r )[0]
    return M_c

def Soliton_r_c_M_c_from_r_c_NumInt( r_c ): # in Msun*kpc
    return r_c*Soliton_M_c_from_r_c_NumInt( r_c, r_c )


Soliton_M_c_from_CHrelation        = 0.25*(M_h/Mmin_0)**(1.0/3.0)*Mmin_0
Soliton_r_c_from_CHrelation        = 1.6/m_22*(M_h/1.0e9)**(-1.0/3.0)
Soliton_rho_c_from_CHrelation      = Soliton_rho_c_from_r_c( Soliton_r_c_from_CHrelation )
Soliton_r_c_M_c_from_CHrelation    = Soliton_r_c_from_CHrelation*Soliton_M_c_from_CHrelation
Soliton_r_c_M_c_from_CHrelation_2  = 400.0/(m_22**2)*(4.4e7)**(2.0/3.0)

r_c_test = Soliton_r_c_from_CHrelation

print("Soliton_rcMc_CHR_1     = %21.14e"%( Soliton_r_c_M_c_from_CHrelation                  ) )
print("Soliton_rcMc_CHR_2     = %21.14e"%( Soliton_r_c_M_c_from_CHrelation_2                ) )
print("Soliton_rcMc_anlytical = %21.14e"%( Soliton_r_c_M_c_from_r_c(        r_c_test )      ) )
print("Soliton_rcMc_NumInt    = %21.14e"%( Soliton_r_c_M_c_from_r_c_NumInt( r_c_test )      ) )
print("Soliton_rcMc_Chiang    = %21.14e"%( Soliton_r_c_M_c_from_r_c_Chiang( r_c_test )      ) )
print("Soliton_rcMs/4_Chiang  = %21.14e"%( Soliton_r_c_M_s_from_r_c_Chiang( r_c_test )*0.25 ) )
print("Soliton_rcMs/4_Chan    = %21.14e"%( Soliton_r_c_M_s_from_r_c_Chan(   r_c_test )*0.25 ) )
print("Soliton_rcMs/4_Chan_2  = %21.14e"%( Soliton_r_c_M_s_from_r_c_Chan_2( r_c_test )*0.25 ) )

print("Soliton_rho_c          = %21.14e"%( Soliton_rho_c_from_r_c(          r_c_test )      ) )
print("Soliton_rho_c_Chan     = %21.14e"%( Soliton_rho_c_from_r_c_Chan(     r_c_test )      ) )
