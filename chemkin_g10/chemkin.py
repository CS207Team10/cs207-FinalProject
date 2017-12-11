import numpy as np
import xml.etree.ElementTree as ET
import chemkin_g10.computation as cp
from chemkin_g10.db import DatabaseOps as dbops
from scipy.integrate import odeint
import matplotlib.pyplot as plt

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



class Simulator:
    """This class represents a simulator for a system of reversible reactions.
    """

    def __init__(self, rsystem, maxTime, numSample=100, timeScale=1e9, eqThreshold=1e-05):
        self.rsystem = rsystem
        self.maxTime = maxTime
        self.numSample = numSample
        self.timeScale = timeScale
        self.eqThreshold = eqThreshold        

    def solveODE(self):
        """Solve the ODE
        INPUTS:
        =======
        t:  float
            total time of simulation
               
        """       
        def fun(concs, t):
            nu = self.rsystem.nu_prod - self.rsystem.nu_react
            rj = cp.progress_rate(self.rsystem.nu_react, self.rsystem.nu_prod, self.rsystem.k, concs, self.rsystem.T, self.rsystem.a, 
                                  self.rsystem.reversibleFlagList, solvingODE=True)
            return np.dot(nu, rj)

        tout = np.linspace(0, self.maxTime/self.timeScale, self.numSample)

        try:
            self.yout = odeint(fun, self.rsystem.concs, tout)
            self.tout = tout
        except ValueError:
            print("ODE solver aborted!")
            raise
        if len(self.yout) != self.numSample:
            raise ValueError("Invalid yout!")

        # get equilibrium point for each time point
        eq_point = [-1 for i in range(len(self.rsystem))]
        eq_diff = [[0 for j in range(len(self.rsystem))] for i in range(self.numSample)] 
        eq_constant = self.rsystem.equilibrium_constant

        for i, concs in enumerate(self.yout):
            if i == 0: continue # there's no product at the beginning
            current_concs = concs.reshape(len(self.rsystem.species), 1)
            reaction_quotient = np.product(current_concs ** self.rsystem.nu_prod, axis=0) / np.product(current_concs ** self.rsystem.nu_react, axis=0)
            
            for j, rq in enumerate(reaction_quotient):
                if not self.rsystem.reversibleFlagList[j]:
                    continue

                eq_diff[i][j] = abs(rq - eq_constant[j]) / eq_constant[j]

                if eq_point[j] != -1:
                    continue

                if eq_diff[i][j] < self.eqThreshold:
                    eq_point[j] = self.tout[i]
                    # eq_point[j] = i

        self.eq_point = eq_point
        self.eq_diff = eq_diff
        return

    def check_equilibrium(self, index, t):
        """Check if the reaction system has reached equilibrium, by comparing
           reaction quotient to reaction coefficient

        INPUTS:
        =======
        index: an integer
               the index of reaction within the reaction system

        t:     an integer
               a time stamp
        
        RETURN:
        =======
        eq: boolean
            if the reaction system has reached equilibrium
        
        """
        if len(self.yout) != self.numSample:
            raise ValueError("Invalid yout!")
        if self.eq_point[index] == -1:
            return False
        return t >= self.eq_point[index]

    def equilibrium_graph(self):
        """Another way to check if the reaction system has reached equilibrium 
 
        RETURN:
        =======
        eq: boolean
            if the reaction system has reached equilibrium
        
        """ 
        slope_diff = (self.yout[-1] - self.yout[-2])/(self.tout[-1]/len(self.tout))

        critical_slope = max(self.yout[-1])/(self.tout[-1])*1e-07
        return all(s < critical_slope for s in slope_diff)



    def plot_specie_all(self):
        """Plot concentration for all species in the reaction system
       
        """ 
        if len(self.yout) != self.numSample:
            raise ValueError("Invalid yout!")
        plt.plot(self.tout, self.yout)
        plt.legend(self.rsystem.species)
        plt.show(block=False)

    def plot_specie(self, index):
        """Plot concentration for one specie in the reaction system

        INPUTS:
        =======
        index: an integer
               the index of specie 
       
        """ 
        if len(self.yout) != self.numSample:
            raise ValueError("Invalid yout!")
        out = np.transpose(self.yout)[index]
        plt.plot(self.tout, out, label = self.rsystem.species[index])
        plt.legend()
        plt.show(block=False)

    def plot_reaction_all(self):
        """Plot (reaction quotient - equilibrium constant) / equilibrium constant 
        for all reactions in the reaction system, to check when each reaction 
        reaches equilibrium
 
        """
        if len(self.yout) != self.numSample:
            raise ValueError("Invalid yout!")
        plt.plot(self.tout[1:], np.sqrt(self.eq_diff[1:]))
        plt.legend([r.reactMeta['id'] for r in self.rsystem.reactionList])
        plt.show(block=False)

# if __name__ == '__main__':
#     T = 900
#     R = 8.314
#     concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
#     rsystem = ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
#     rsystem.buildFromXml("../tests/data/xml/rxns_reversible.xml", concs)
#     # print("Progress rate: \n", rsystem.getProgressRate(), "\n")
#     print("Reaction rate: \n", rsystem.getReactionRate(), "\n")
#     # print("System info: \n", rsystem, "\n")


#     sim = Simulator(rsystem, 0.05)
#     sim.solveODE()
#     # print(sim.yout)
#     # sim.plot_specie_all()
#     # print(sim.check_equilibrium(5, 5e-11))
#     # sim.plot_specie(4)
#     # print(sim.eq_diff)
#     # sim.plot_reaction_all()
#     # print(sim.equilibrium_graph())










