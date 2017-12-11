import numpy as np
from chemkin_g10.thermo import backward_coeffs

import os
path = os.path.dirname(os.path.realpath(__file__)) + "/../tests/"

def k_const(k=1.0):
    """Simply returns a constant reaction rate coefficient

    INPUTS:
    =======
    k: float, default value = 1.0
       Constant reaction rate coefficient

    RETURNS:
    ========
    k: float
       Constant reaction rate coefficient

    EXAMPLES:
    =========
    >>> k_const(5.0)
    5.0
    """
    if k < 0:
        raise ValueError("Negative reaction rate coefficients are prohibited.")

    return k


def k_arr(A, E, T, R=8.314):
    """Calculates the Arrhenius reaction rate coefficient

    INPUTS:
    =======
    A: float
       Arrhenius prefactor
       Must be positive
    E: float
       Activation energy
    T: float
       Temperature
       Must be positive
    R: float, default value = 8.314
       Ideal gas constant
       Must be positive

    RETURNS:
    ========
    k: float
       Arrhenius reaction rate coefficient

    EXAMPLES:
    =========
    >>> k_arr(2.0, 3.0, 100.0)
    1.9927962618542914
    """

    if A < 0.0:
        raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

    if T < 0.0:
        raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

    if R < 0.0:
        raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

    return A * np.exp(-E / R / T)

def k_mod_arr(A, b, E, T, R=8.314):
    """Calculates the modified Arrhenius reaction rate coefficient

    INPUTS:
    =======
    A: float
       Arrhenius prefactor
       Must be positive
    b: float
       Modified Arrhenius parameter
    E: float
       Activation energy
    T: float
       Temperature
       Must be positive
    R: float, default value = 8.314
       Ideal gas constant
       Must be positive

    RETURNS:
    ========
    k: float
       Modified Arrhenius reaction rate coefficient

    EXAMPLES:
    =========
    >>> k_mod_arr(2.0, -0.5, 3.0, 100.0)
    0.19927962618542916
    """
    if A < 0.0:
        raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

    if T < 0.0:
        raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

    if R < 0.0:
        raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

    return A * T**b * np.exp(-E / R / T)

def progress_rate(nu_react, nu_prod, k, concs, T, a, reversibleFlagList, solvingODE=False):
    """Returns the progress rate of a system elementary reactions (whether reversible or not)
    INPUTS:
    =======
    nu_react: numpy array of floats,
              size: num_species X num_reactions
              stoichiometric coefficients for the reaction
    nu_prod:  numpy array of floats,
              size: num_species X num_reactions
              stoichiometric coefficients for the products
    k:        array of floats
              Reaction rate coefficient for the reaction
    concs:    numpy array of floats
              concentration of species
    T:        float
              Temperature
              Must be positive
    a:        numpy array of floats,
              size: num_species X 7
              nasa coefficients of each species
    reversibleFlagList:
              array of booleans, 
              size: num_reactions
              boolean indicator indicating wether each reaction is reversible
    RETURNS:
    ========
    progress: numpy array of floats
              size: num_reactions
              progress rate of each reaction
    
    """
    progress = k.copy() # Initialize progress rates with reaction rate coefficients

    for jdx, rj in enumerate(progress):
        if rj < 0:
            raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))

        # initialize the subtraction part for reversible reaction
        if reversibleFlagList[jdx] == True:
            sub = backward_coeffs(rj, nu_prod[:, jdx] - nu_react[:, jdx], T, a)
        else:
            sub = 0

        for idx, xi in enumerate(concs):
            nu_ij_r = nu_react[idx,jdx]
            nu_ij_p = nu_prod[idx,jdx]
            if not solvingODE and xi  < 0.0:
                raise ValueError("x{0} = {1:18.16e}:  Negative concentrations are prohibited!".format(idx, xi))
            if nu_ij_r < 0:
                raise ValueError("nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx, nu_ij_r))
            if nu_ij_p < 0:
                raise ValueError("nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx, nu_ij_p))

            # forward progress rate
            progress[jdx] *= xi**nu_ij_r

            # calculate the backward rate if this reaction is reversible
            if reversibleFlagList[jdx] == True:
                sub *= xi**nu_ij_p

        # subtract it to get the overall progress rate
        progress[jdx] -= sub

    return progress

def reaction_rate(nu_react, nu_prod, k, concs, T, a, reversibleFlagList):
    """Returns the progress rate of a system elementary reactions (whether reversible or not)
    INPUTS:
    =======
    nu_react: numpy array of floats,
              size: num_species X num_reactions
              stoichiometric coefficients for the reaction
    nu_prod:  numpy array of floats,
              size: num_species X num_reactions
              stoichiometric coefficients for the products
    k:        array of floats
              Reaction rate coefficient for the reaction
    concs:    numpy array of floats
              concentration of species
    T:        float
              Temperature
              Must be positive
    a:        numpy array of floats,
              size: num_species X 7
              nasa coefficients of each species
    reversibleFlagList:
              array of booleans, 
              size: num_reactions
              boolean indicator indicating wether each reaction is reversible
    RETURNS:
    ========
    f: numpy array of floats
       size: num_species
       reaction rate of each species
    
    """
    nu = nu_prod - nu_react
    rj = progress_rate(nu_react, nu_prod, k, concs, T, a, reversibleFlagList)
    return np.dot(nu, rj)


def equilibrium_constant(nu_react, nu_prod, k, T, a, reversibleFlagList):
    """Returns the progress rate of a system elementary reactions (whether reversible or not)
    INPUTS:
    =======
    nu_react: numpy array of floats,
              size: num_species X num_reactions
              stoichiometric coefficients for the reaction
    nu_prod:  numpy array of floats,
              size: num_species X num_reactions
              stoichiometric coefficients for the products
    k:        array of floats
              Reaction rate coefficient for the reaction
    T:        float
              Temperature
              Must be positive
    a:        numpy array of floats,
              size: num_species X 7
              nasa coefficients of each species
    reversibleFlagList:
              array of booleans, 
              size: num_reactions
              boolean indicator indicating wether each reaction is reversible
    RETURNS:
    ========
    f: numpy array of floats
       size: num_reactions
       the equilibrium constant of each reaction
    
    """
    eq_constant = []
    for jdx, kj in enumerate(k):
        if kj < 0:
            raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(kj))

        # initialize the subtraction part for reversible reaction
        if reversibleFlagList[jdx] == True:
            kb = backward_coeffs(kj, nu_prod[:, jdx] - nu_react[:, jdx], T, a)
            eq_constant.append(kj / kb) # ke
        else:
            eq_constant.append(0) # no such constant for irreversible reaction

    return eq_constant



if __name__ == '__main__':
    print("")

