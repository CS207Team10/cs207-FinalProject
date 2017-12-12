import numpy as np
import xml.etree.ElementTree as ET
import chemkin_g10.computation as cp
from chemkin_g10.db import DatabaseOps as dbops

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
            self.k = cp.k_const(meta['k'])
        elif self.rateCoeffMeta['type'] =="Arrhenius":
            self.k = cp.k_arr(meta['A'], meta['E'], meta['T'], meta['R'])
        elif self.rateCoeffMeta['type'] =="modifiedArrhenius":
            self.k = cp.k_mod_arr(meta['A'], meta['b'], meta['E'], meta['T'], meta['R'])
        else:
            # Other type of reaction rate coeff
            self.k = None # k = cp.newMethodToComputeK(...)
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
    dbFileName: string
                the name of the db file
    """
    def __init__(self, T, R, dbFileName):
        self.T = T
        self.R = R
        self.dbops = dbops(dbFileName)


    def buildFromList(self, reactionList, species, concs):
        """Build ReactionSystem from a reactionList.

        INPUTS:
        =======
        reactionList: list of Reaction object
        species:      list of str
        concs:        list of float
                      initial concerntration for each species

        """
        if len(species) != len(concs):
            raise ValueError("Size of concentration does not match to number of species!")

        self.concs = concs
        self.reactionList = reactionList
        self.nu_react = np.array([r.reactCoeff for r in self.reactionList]).T
        self.nu_prod = np.array([r.productCoeff for r in self.reactionList]).T
        self.k = np.array([r.k for r in self.reactionList])
        self.a = self.dbops.get_coeffs(self.species, self.T)
        self.reversibleFlagList = [r.reactMeta['reversible']=='yes' for r in reactionList]
        self.progress_rate = cp.progress_rate(self.nu_react, self.nu_prod, self.k, self.concs, self.T, self.a, self.reversibleFlagList)
        self.reaction_rate = cp.reaction_rate(self.nu_react, self.nu_prod, self.k, self.concs, self.T, self.a, self.reversibleFlagList)
        self.equilibrium_constant = cp.equilibrium_constant(self.nu_react, self.nu_prod, self.k, self.T, self.a, self.reversibleFlagList)


    def buildFromXml(self, inputFile, concs):
        """Build ReactionSystem from an input file.

        INPUTS:
        =======
        inputFile: str
                   file name of an XML file

        """

        reactionList, species = self.parse(inputFile, self.T, self.R)
        if len(species) != len(concs):
            raise ValueError("Size of concentration does not match to number of species!")

        self.inputFile = inputFile
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
        >>> type( cp.parse( "./test1.xml", 340, 8.314)[0] )
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
                k = cp.k_const(float(coeffSection.find("Constant").find("k").text))
            elif coeffSection.find("Arrhenius") != None:
                rateCoeffMeta["type"] = "Arrhenius"
                rateCoeffMeta["A"] = float(coeffSection.find("Arrhenius").find("A").text)
                rateCoeffMeta["E"] = float(coeffSection.find("Arrhenius").find("E").text)
                k = cp.k_arr(rateCoeffMeta["A"], rateCoeffMeta["E"], T, R)
            elif coeffSection.find("modifiedArrhenius") != None:
                rateCoeffMeta["type"] = "modifiedArrhenius"
                rateCoeffMeta["A"] = float(coeffSection.find("modifiedArrhenius").find("A").text)
                rateCoeffMeta["b"] = float(coeffSection.find("modifiedArrhenius").find("b").text)
                rateCoeffMeta["E"] = float(coeffSection.find("modifiedArrhenius").find("E").text)
                k = cp.k_mod_arr(rateCoeffMeta["A"], rateCoeffMeta["b"], rateCoeffMeta["E"], T, R)
            else:
                # Other type of reaction rate coeff
                k = None # k = cp.newMethodToComputeK(...)

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

    def __str__(self):
        res = "\n"
        res += "The system:"
        for r in self.reactionList:
            res += "\n" + str(r)
        return res + "\n"

    def __len__(self):
        return len(self.reactionList)


# if __name__ == '__main__':
#     T = 900
#     R = 8.314
#     concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
#     rsystem = ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
#     rsystem.buildFromXml("../tests/data/xml/rxns_reversible.xml", concs)
#     print("Reaction rate: \n", rsystem.getReactionRate(), "\n")
#     # print("Progress rate: \n", rsystem.getProgressRate(), "\n")
#     # print("System info: \n", rsystem, "\n")
