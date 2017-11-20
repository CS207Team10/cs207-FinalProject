import numpy as np

def H_over_RT(T, a):

    # WARNING:  This line will depend on your own data structures!
    # Be careful to get the correct coefficients for the appropriate 
    # temperature range.  That is, for T <= Tmid get the low temperature 
    # range coeffs and for T > Tmid get the high temperature range coeffs.
    H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
            + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
            + a[:,5] / T)

    return H_RT
           
def S_over_R(T, a):

    # WARNING:  This line will depend on your own data structures!
    # Be careful to get the correct coefficients for the appropriate 
    # temperature range.  That is, for T <= Tmid get the low temperature 
    # range coeffs and for T > Tmid get the high temperature range coeffs.
    S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
           + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

    return S_R

def backward_coeffs(kf, nu, T, a, p0=100000, R=8.3144598):

    # Change in enthalpy and entropy for each reaction
    delta_H_over_RT = np.dot(nu.T, H_over_RT(T, a))
    delta_S_over_R = np.dot(nu.T, S_over_R(T, a))

    # Negative of change in Gibbs free energy for each reaction 
    delta_G_over_RT = delta_S_over_R - delta_H_over_RT

    # Prefactor in Ke
    fact = p0 / R / T

    # Ke
    gamma = np.sum(nu, axis=0)
    kb = fact**gamma * np.exp(delta_G_over_RT)

    return kf / kb

