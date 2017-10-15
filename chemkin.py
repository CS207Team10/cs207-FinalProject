import numpy as np
import xml.etree.ElementTree as ET

# Some globals
T_DEFAULT = 1500
R_DEFAULT = 8.314

class ChemUtil:

    @classmethod
    def parse(cls, inputFile, T, R):
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
            # print(reactStr)

            # Parse reaction rate coeff parameters
            rateCoeffType = None
            k = None
            rateCoeffPar = []
            coeffSection = row.find("rateCoeff")
            if coeffSection.find("Constant") != None:
                rateCoeffType = "Constant"
                rateCoeffPar.append(float(coeffSection.find("Constant").find("k").text))
                k = cls.k_const(rateCoeffPar[0])
            elif coeffSection.find("Arrhenius") != None:
                rateCoeffType = "Arrhenius"
                rateCoeffPar.append(float(coeffSection.find("Arrhenius").find("A").text))
                rateCoeffPar.append(float(coeffSection.find("Arrhenius").find("E").text))
                k = cls.k_arr(rateCoeffPar[0], rateCoeffPar[1], T, R)
            else:
                rateCoeffType = "modifiedArrhenius"
                rateCoeffPar.append(float(coeffSection.find("modifiedArrhenius").find("A").text))
                rateCoeffPar.append(float(coeffSection.find("modifiedArrhenius").find("b").text))
                rateCoeffPar.append(float(coeffSection.find("modifiedArrhenius").find("E").text))
                k = cls.k_mod_arr(rateCoeffPar[0], rateCoeffPar[1], rateCoeffPar[2], T, R)

            # Put coeffs of reactants, products
            reactCoeff = np.zeros(len(species))
            productCoeff = np.zeros(len(species))
            # Split the "_:_" pairs and get the coeff for corresponding species
            for mp in row.find("reactants").text.strip().split(" "):
                sp, co = mp.split(":")
                reactCoeff[species.index(sp)] = float(co)
            for mp in row.find("products").text.strip().split(" "):
                sp, co = mp.split(":")
                productCoeff[species.index(sp)] = float(co)

            # save the metadata (reversible, type, id ...)
            reactMeta = row.attrib
            reactMeta["T"] = T
            reactMeta["R"] = R
            # Create the Reaction object, and add to the list

            reaction = Reaction(reactStr, k, rateCoeffType, rateCoeffPar,
                                reactCoeff, productCoeff, reactMeta)
            reactionList.append(reaction)
            # print(reactStr, k, rateCoeffType, rateCoeffPar, reactCoeff, productCoeff, reactMeta)

        return reactionList

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
        >>> k_const(5.0)
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

    @classmethod
    def progress_rate(cls, nu_react, k, concs):
        """Returns the progress rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        nu_react: numpy array of floats,
                  size: num_species X num_reactions
                  stoichiometric coefficients for the reaction
        k:        array of floats
                  Reaction rate coefficient for the reaction
        concs:    numpy array of floats
                  concentration of species

        RETURNS:
        ========
        omega: numpy array of floats
               size: num_reactions
               progress rate of each reaction

        EXAMPLES:
        =========
        >>> progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([2.0, 1.0, 1.0]), 10.0)
        array([ 40.,  20.])
        """
        progress = k.copy() # Initialize progress rates with reaction rate coefficients
        for jdx, rj in enumerate(progress):
            if rj < 0:
                raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))
            for idx, xi in enumerate(concs):
                nu_ij = nu_react[idx,jdx]
                if xi  < 0.0:
                    raise ValueError("x{0} = {1:18.16e}:  Negative concentrations are prohibited!".format(idx, xi))
                if nu_ij < 0:
                    raise ValueError("nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx, nu_ij))

                progress[jdx] *= xi**nu_ij
        return progress

    @classmethod
    def reaction_rate(cls, nu_react, nu_prod, k, concs):
        """Returns the reaction rate of a system of irreversible, elementary reactions

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
        ...to be added
        """
        nu = nu_prod - nu_react
        rj = cls.progress_rate(nu_react, k, concs)
        return np.dot(nu, rj)


class Reaction:

    def __init__(self, reactStr,
                 k, rateCoeffType, rateCoeffPar,
                 reactCoeff, productCoeff,
                 reactMeta):
        self.reactStr = reactStr
        self.k = k
        self.rateCoeffType = rateCoeffType
        self.rateCoeffPar = rateCoeffPar
        self.reactCoeff = reactCoeff
        self.productCoeff = productCoeff
        self.reactMeta = reactMeta
        # self.setRateCoeff()

    # def setMeta(self, **args): # used to set T, R, reversible, elementary, id
    #     for par in args:
    #         self.reactMeta[par] = args[par]
    #
    # def setRateCoeff(self): # set the rate coeff
    #     k = None
    #     T = self.reactMeta["T"]
    #     R = self.reactMeta["R"]
    #     par = self.rateCoeffPar
    #     if self.rateCoeffType == "Constant":
    #         k = self.k_const(par[0])
    #     elif self.rateCoeffType == "Arrhenius":
    #         k = self.k_arr(par[0], par[1], T, R)
    #     elif self.rateCoeffType == "modifiedArrhenius":
    #         k = self.k_arr(par[0], par[1], par[2], T, R)
    #     # others would cause errors during parsing
    #
    #     self.k = k


    def __str__(self):
        return self.reactStr




class ReactionSystem:
    def __init__(self, inputFile, T, R, concs):
        self.inputFile = inputFile
        self.T = T
        self.R = R
        self.concs = concs

    def build(self):
        self.reactionList = ChemUtil.parse(self.inputFile, self.T, self.R)
        self.nu_react = np.array([r.reactCoeff for r in self.reactionList]).T
        self.nu_prod = np.array([r.productCoeff for r in self.reactionList]).T
        self.k = np.array([r.k for r in self.reactionList])
        # print(self.nu_react, self.nu_prod)
        self.progress_rate = ChemUtil.progress_rate(self.nu_react, self.k, self.concs)
        self.reaction_rate = ChemUtil.reaction_rate(self.nu_react, self.nu_prod, self.k, self.concs)



if __name__ == '__main__':
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ReactionSystem("test1.xml", T_DEFAULT, R_DEFAULT, concs)
    rsystem.build()

    print(rsystem.progress_rate)
    print(rsystem.reaction_rate)
