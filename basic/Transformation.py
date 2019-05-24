import numpy as np

import MKastro.basic.Constants as cst


#################################
##    Flux Transformations     ##
#################################

def T_to_R( T, Z, A, unit_m='GeV' ):
    m = A * cst.m_proton[unit_m]
    E = T + m
    p = np.sqrt(E*E-m*m)
    return p/np.fabs(1.*Z)

def R_to_T( R, Z, A, unit_m='GeV' ):
    m = A * cst.m_proton[unit_m]
    p = R*np.fabs(1.*Z)
    return np.sqrt(p*p+m*m)-m

def Tn_to_R( Tn, Z, A, unit_m='GeV' ):
    m_p = cst.m_proton[unit_m]
    #E   = Tn*A + m
    p   = np.sqrt(Tn*Tn+2*Tn*m_p)*A
    return p/np.fabs(1.*Z)

def R_to_Tn( R, Z, A, unit_m='GeV' ):
    m = A * cst.m_proton[unit_m]
    p = R*np.fabs(1.*Z)
    return (np.sqrt(p*p+m*m)-m)/A

def dT_by_dR( T, Z, A, unit_m='GeV' ):
    m = A * cst.m_proton[unit_m]
    E = T + m
    p = np.sqrt(E*E-m*m)
    return p*np.fabs(1.*Z)/E

def dR_by_dT( R,  Z, A, unit_m='GeV'):
    m = A * cst.m_proton[unit_m]
    p = R*np.fabs(1.*Z)
    E = np.sqrt(p*p+m*m)
    return E/(p*np.fabs(1.*Z))

def ConvertSpectrumFromTnToR( Tn, F, Z, A, unit_m='GeV' ):
    R = T_to_R(Tn*A, Z, A)
    Fnew = F*dT_by_dR(Tn*A,Z,A)/A;
    return R, Fnew

def ConvertSpectrumFromRToTn( R, F, Z, A, unit_m='GeV' ):
    Tn = R_to_T(R, Z, A)/A;
    Fnew = F*dR_by_dT(R, Z, A)*A;
    return Tn, Fnew

def SolarModulationTn( Tn, F, phi, Z, A, unit_m='GeV'):
    Tn_new    = Tn - np.fabs(1.*Z)/A*phi;
    F_new           = F*(Tn_new*Tn_new+2*Tn_new*cst.m_proton[unit_m])/(Tn*Tn+2*Tn*cst.m_proton[unit_m])
    return Tn_new, F_new

def SolarModulationTn_interpolated( Tn_res, Tn, F, phi, Z, A, unit_m='GeV'):
    Tn_original = Tn_res + np.fabs(1.*Z)/A*phi
    F_intp      = np.exp( np.interp( np.log(Tn_original), np.log(Tn), np.log(F+1e-44) ) )
    F_new           = F_intp*(Tn_res**2+2*Tn_res*cst.m_proton[unit_m])/(Tn_original**2+2*Tn_original*cst.m_proton[unit_m])
    return  F_new

def SolarModulationR( R, Fr, phi, Z, A, unit_m='GeV'):
    Tn, F = ConvertSpectrumFromRToTn(R, Fr, Z, A)
    Tn_new    = Tn - np.fabs(1.*Z)/A*phi;
    F_new           = F*(Tn_new*Tn_new+2*Tn_new*cst.m_proton[unit_m])/(Tn*Tn+2*Tn*cst.m_proton[unit_m])
    R_new, Fr_new = ConvertSpectrumFromTnToR(Tn_new, F_new, Z, A)
    return R_new, Fr_new
