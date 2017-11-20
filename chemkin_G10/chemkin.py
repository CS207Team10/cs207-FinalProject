import numpy as np
import xml.etree.ElementTree as ET
import sqlite3
import os
from pathlib import Path

# Some globals for reference
T_DEFAULT = 1500
R_DEFAULT = 8.314
path = os.path.dirname(os.path.realpath(__file__)) + "/../tests/"

class ChemUtil:
    """The class that contains all the utility functions"""

    @classmethod
    def k_const(cls, k=1.0):
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
        >>> ChemUtil.k_const(5.0)
        5.0
        """
        if k < 0:
            raise ValueError("Negative reaction rate coefficients are prohibited.")

        return k

    @classmethod
    def k_arr(cls, A, E, T, R=8.314):
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
        >>> ChemUtil.k_arr(2.0, 3.0, 100.0)
        1.9927962618542914
        """

        if A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

        return A * np.exp(-E / R / T)

    @classmethod
    def k_mod_arr(cls, A, b, E, T, R=8.314):
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
        >>> ChemUtil.k_mod_arr(2.0, -0.5, 3.0, 100.0)
        0.19927962618542916
        """
        if A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(R))

        return A * T**b * np.exp(-E / R / T)

    @classmethod
    def progress_rate(cls, nu_react, nu_prod, k, concs, T, a, reversibleFlagList):
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

        RETURNS:
        ========
        progress: numpy array of floats
                  size: num_reactions
                  progress rate of each reaction

        EXAMPLES:
        =========
        >>> ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]))
        array([ 40.,  20.])
        """
        progress = k.copy() # Initialize progress rates with reaction rate coefficients

        for jdx, rj in enumerate(progress):
            if rj < 0:
                raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))

            # initialize the subtraction part for reversible reaction
            if reversibleFlagList[jdx] == True:
                sub = ChemUtil.backward_coeffs(rj, nu_prod[:, jdx] - nu_react[:, jdx], T, a)
            else:
                sub = 0

            for idx, xi in enumerate(concs):
                nu_ij_r = nu_react[idx,jdx]
                nu_ij_p = nu_prod[idx,jdx]
                if xi  < 0.0:
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

    @classmethod
    def reaction_rate(cls, nu_react, nu_prod, k, concs, T, a, reversibleFlagList):
        """Returns the reaction rate of a system of elementary reactions

        INPUTS:
        =======
        nu_react: numpy array of floats,
                  size: num_species X num_reactions
                  stoichiometric coefficients for the reactants
        nu_prod:  numpy array of floats,
                  size: num_species X num_reactions
                  stoichiometric coefficients for the products
        k:        float, default value is 10,
                  Reaction rate coefficient for the reaction
        concs:    numpy array of floats
                  concentration of species

        RETURNS:
        ========
        f: numpy array of floats
           size: num_species
           reaction rate of each specie

        EXAMPLES:
        =========
        >>> ChemUtil.reaction_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([10.0, 10.0]), np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0])) 
        array([ 500.,  560.,  580.])
        """
        nu = nu_prod - nu_react
        rj = cls.progress_rate(nu_react, nu_prod, k, concs, T, a, reversibleFlagList)
        return np.dot(nu, rj)

    @classmethod
    def H_over_RT(cls, T, a):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
                + a[:,5] / T)

        return H_RT
               
    @classmethod
    def S_over_R(cls, T, a):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

        return S_R

    @classmethod
    def backward_coeffs(cls, kf, nu, T, a, p0=100000, R=8.3144598):

        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(nu.T, cls.H_over_RT(T, a))
        delta_S_over_R = np.dot(nu.T, cls.S_over_R(T, a))

        # Negative of change in Gibbs free energy for each reaction 
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = p0 / R / T

        # Ke
        gamma = np.sum(nu, axis=0)
        kb = fact**gamma * np.exp(delta_G_over_RT)

        return kf / kb


    @classmethod
    def get_nasa_coeffs(cls, cursor, species, T):
        a = []
        for s in species:
            # check table LOW

            query = '''SELECT * FROM LOW WHERE SPECIES_NAME="{}"'''.format(s)
            res = cursor.execute(query).fetchall()
            if len(res) != 1:
                raise ValueError("Specie {} not in the database!".format(s))
            if res[0][2] <= T and  T <= res[0][3]: # in the range
                a.append(list(res[0][-7:]))
                continue

            # check table HIGH if not in the LOW's range
            query = '''SELECT * FROM HIGH WHERE SPECIES_NAME="{}"'''.format(s)
            res = cursor.execute(query).fetchall()
            if len(res) != 1:
                raise ValueError("Specie {} not in the database!".format(s))
            if res[0][2] <= T and  T <= res[0][3]:
                a.append(list(res[0][-7:]))
                continue

            # The temperature T is beyond the range
            raise ValueError("The temperate is beyond the range for species {}!".format(s))
        return np.array(a)


    @classmethod
    def parse(cls, inputFile, T, R):
        """Parse an XML input file and return list of chemkin.Reaction metadata

        INPUTS:
        =======
        inputFile:XML file,
        T:        float
                  Temperature
                  Must be positive
        R:        float,
                  Ideal gas constant
                  Must be positive

        RETURNS:
        ========
        reactionList: list of chemkin.Reaction metadata
                      
        EXAMPLES:
        =========
        >>> type( ChemUtil.parse( "./test1.xml", 340, 8.314)[0] )   
        <class 'chemkin.Reaction'>
        """
        try:
            tree = ET.parse(inputFile)
        except:
            raise ValueError("Must be a valid xml file!")
        root = tree.getroot()
        species = root.find("phase").find("speciesArray").text.strip().split(" ")
        # print(species)

        reactionData = root.find("reactionData")
        reactionList = [] # list of "Reaction"
        for row in reactionData:

            # reaction formula
            reactStr = row.find("equation").text

            # metadata for rate coeff and reaction
            rateCoeffMeta = dict()
            rateCoeffMeta["T"] = T
            rateCoeffMeta["R"] = R
            reactMeta = row.attrib # reversible/irreversible, type, id ...

            # Parse reaction rate coeff parameters and save to rateCoeffMeta
            coeffSection = row.find("rateCoeff")
            k = None
            if coeffSection.find("Constant") != None:
                rateCoeffMeta["type"] = "Constant"
                k = cls.k_const(float(coeffSection.find("Constant").find("k").text))
            elif coeffSection.find("Arrhenius") != None:
                rateCoeffMeta["type"] = "Arrhenius"
                rateCoeffMeta["A"] = float(coeffSection.find("Arrhenius").find("A").text)
                rateCoeffMeta["E"] = float(coeffSection.find("Arrhenius").find("E").text)
                k = cls.k_arr(rateCoeffMeta["A"], rateCoeffMeta["E"], T, R)
            elif coeffSection.find("modifiedArrhenius") != None:
                rateCoeffMeta["type"] = "modifiedArrhenius"
                rateCoeffMeta["A"] = float(coeffSection.find("modifiedArrhenius").find("A").text)
                rateCoeffMeta["b"] = float(coeffSection.find("modifiedArrhenius").find("b").text)
                rateCoeffMeta["E"] = float(coeffSection.find("modifiedArrhenius").find("E").text)
                k = cls.k_mod_arr(rateCoeffMeta["A"], rateCoeffMeta["b"], rateCoeffMeta["E"], T, R)
            else:
                # Other type of reaction rate coeff
                k = None # k = ChemUtil.newMethodToComputeK(...)

            # Coeffs of reactants, products
            # Split the "_:_" pairs and get the coeff for corresponding species
            reactCoeff = np.zeros(len(species))
            productCoeff = np.zeros(len(species))
            for mp in row.find("reactants").text.strip().split(" "):
                sp, co = mp.split(":")
                reactCoeff[species.index(sp)] = float(co)
            for mp in row.find("products").text.strip().split(" "):
                sp, co = mp.split(":")
                productCoeff[species.index(sp)] = float(co)

            # Create the Reaction object, and add to the list
            reaction = Reaction(reactStr, k, reactCoeff, productCoeff,
                                rateCoeffMeta, reactMeta)
            reactionList.append(reaction)

        return reactionList, species


class Reaction:
    """The class that represent a single Reaction

    Parameters
    ----------
    reactStr:       str
                    string representation of the reaction
    k:              float, default value is 10,
                    Reaction rate coefficient for the reaction
    reactCoeff:     numpy array of floats,
                    size: num_species X num_reactions
                    stoichiometric coefficients for the reactants
    productCoeff:   numpy array of floats,
                    size: num_species X num_reactions
                    stoichiometric coefficients for the products
    rateCoeffMeta:  dict
                    metadata for the reaction rate coefficient
    reactMeta:      dict
                    metadata for the reaction
    """
    def __init__(self, reactStr,
                 k, reactCoeff, productCoeff,
                 rateCoeffMeta, reactMeta):
        self.reactStr = reactStr
        self.k = k
        self.reactCoeff = reactCoeff
        self.productCoeff = productCoeff
        self.rateCoeffMeta = rateCoeffMeta
        self.reactMeta = reactMeta

    def updateCoeff(self, **args):
        """update the metadata of reaciton rate coefficient and 
           recalculate the coefficient.
           
        INPUTS:
        =======
        args: T=..., R=..., type=..., A=..., b=..., E=... 
        """
        for par in args:
            self.rateCoeffMeta[par] = args[par]
        meta = self.rateCoeffMeta
        if self.rateCoeffMeta['type'] =="constant":
            self.k = ChemUtil.k_const(meta['k'])
        elif self.rateCoeffMeta['type'] =="Arrhenius":
            self.k = ChemUtil.k_arr(meta['A'], meta['E'], meta['T'], meta['R'])
        elif self.rateCoeffMeta['type'] =="modifiedArrhenius":
            self.k = ChemUtil.k_mod_arr(meta['A'], meta['b'], meta['E'], meta['T'], meta['R'])
        else:
            # Other type of reaction rate coeff
            self.k = None # k = ChemUtil.newMethodToComputeK(...)
        return

    def updateReaction(self, **args):
        """Update the metadata of the reactionthe possible parameters: "duplicate reactions"
           
        INPUTS:
        =======
        args: reversible=..., type=..., id=...
        """
        for par in args:
            self.reactMeta[par] = args[par]
        # Future use, do something
        return

    def __str__(self):
        return self.reactStr + "\n" + \
               "\t" + "rate coeff: " + str(self.k) + "\n" + \
               "\t" + "rate coeff metadata: " + str(self.rateCoeffMeta) + "\n" + \
               "\t" + "reaction metadata: " + str(self.reactMeta)


class ReactionSystem:
    """The class that represents a system of reactions

    Parameters
    ----------
    T:     float or int
           the temperature of the system
    R:     float or int
           the universal gas constant
    concs: numpy array of floats
           concentration of species
    """
    def __init__(self, T, R, dbFileName):
        self.T = T
        self.R = R
        my_file = Path(path + "/data/db/" + dbFileName)
        if my_file.is_file():
            self.db = sqlite3.connect(str(my_file))
        else:
            raise ValueError("The db file does not exist!")


    def buildFromList(self, reactionList, species, concs):
        """Build ReactionSystem from a reactionList.

        INPUTS:
        =======
        reactionList: list of Reaction object
        
        """
        if len(species) != len(concs):
            raise ValueError("Size of concentration does not match to number of species!")

        self.concs = concs
        self.reactionList = reactionList
        self.nu_react = np.array([r.reactCoeff for r in self.reactionList]).T
        self.nu_prod = np.array([r.productCoeff for r in self.reactionList]).T
        self.k = np.array([r.k for r in self.reactionList])
        self.a = ChemUtil.get_nasa_coeffs(self.db.cursor(), self.species, self.T)
        reversibleFlagList = [r.reactMeta['reversible']=='yes' for r in reactionList]
        self.progress_rate = ChemUtil.progress_rate(self.nu_react, self.nu_prod, self.k, self.concs, self.T, self.a, reversibleFlagList)
        self.reaction_rate = ChemUtil.reaction_rate(self.nu_react, self.nu_prod, self.k, self.concs, self.T, self.a, reversibleFlagList)


    def buildFromXml(self, inputFile, concs):
        """Build ReactionSystem from an input file.

        INPUTS:
        =======
        inputFile: str
                   file name of an XML file
        
        """
        
        reactionList, species = ChemUtil.parse(inputFile, self.T, self.R)
        if len(species) != len(concs):
            raise ValueError("Size of concentration does not match to number of species!")

        self.reactionList = reactionList
        self.species = species
        self.concs = concs
        self.buildFromList(self.reactionList, self.species, self.concs)

    def getProgressRate(self):
        """Return progress rate.

        RETURN:
        =======
        progress_rate: progress rate of the reaction system
        
        """
        return self.progress_rate

    def getReactionRate(self):
        """Return reaction rate.

        RETURN:
        =======
        reaction_rate: reaction rate of the reaction system
        
        """
        return self.reaction_rate

    def __str__(self):
        res = "\n"
        res += "The system:"
        for r in self.reactionList:
            res += "\n" + str(r)
        return res + "\n"

    def __len__(self):
        return len(self.reactionList)


if __name__ == '__main__':

    concs = np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    rsystem = ReactionSystem(T_DEFAULT, R_DEFAULT, "nasa.sqlite")
    rsystem.buildFromXml(path+ "/data/xml/rxns_reversible.xml", concs)
    print(rsystem.getProgressRate())
    print(rsystem.getReactionRate())




