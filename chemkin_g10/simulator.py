import numpy as np
import chemkin_g10.computation as cp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import requests
import json
import webbrowser

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
        if not hasattr(self, 'yout'):
            raise ValueError("Please solve ODE first!")
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
        if not hasattr(self, 'yout'):
            raise ValueError("Please solve ODE first!")
        slope_diff = (self.yout[-1] - self.yout[-2])/(self.tout[-1]/len(self.tout))

        critical_slope = max(self.yout[-1])/(self.tout[-1])*1e-07
        return all(s < critical_slope for s in slope_diff)


    def plot_specie_all(self):
        """Plot concentration for all species in the reaction system

        """
        if not hasattr(self, 'yout'):
            raise ValueError("Please solve ODE first!")
        plt.plot(self.tout, self.yout)
        plt.legend(self.rsystem.species)
        plt.show()

    def plot_specie(self, species):
        """Plot concentration for one specie in the reaction system

        INPUTS:
        =======
        species: string
                 certain species
        """
        if not hasattr(self, 'yout'):
            raise ValueError("Please solve ODE first!")
        if species not in self.rsystem.species:
            raise ValueError("Please provide a valid species!")
        index = self.rsystem.species.index(species)
        out = np.transpose(self.yout)[index]
        plt.plot(self.tout, out, label = self.rsystem.species[index])
        plt.legend()
        plt.show()

    def plot_reaction_all(self):
        """Plot (reaction quotient - equilibrium constant) / equilibrium constant
        for all reactions in the reaction system, to check when each reaction
        reaches equilibrium

        """
        if not hasattr(self, 'yout'):
            raise ValueError("Please solve ODE first!")
        plt.plot(self.tout[1:], np.sqrt(self.eq_diff[1:]))
        plt.legend([r.reactMeta['id'] for r in self.rsystem.reactionList])
        plt.show()

    # def visualize_dev(self):
    #     """Return a link to visualize the system in the web front-end (javascript), which is
    #        way better than tkinter
    #     """
    #     url_local = "http://127.0.0.1:8000/visualization/viz/"
    #     files = {'file': open(self.rsystem.inputFile, 'rb')}
    #     r = requests.post(url_local, data = {'T':self.rsystem.T, 'concs':json.dumps(self.rsystem.concs.tolist()), 'maxTime':self.maxTime}, files=files)
    #     url_viz = "http://127.0.0.1:8000/visualization/demo/" + r.text + "/"
    #     webbrowser.open_new(url_viz)
    #     return

    def visualize(self):
        """Return a link to visualize the system in the web front-end (javascript) instead of using tkinter, as it provides better visualizatons.
        """
        url_local = "http://cs207g10viz.us-east-1.elasticbeanstalk.com/visualization/viz/"
        files = {'file': open(self.rsystem.inputFile, 'rb')}
        r = requests.post(url_local, data = {'T':self.rsystem.T, 'concs':json.dumps(self.rsystem.concs.tolist()), 'maxTime':self.maxTime}, files=files)
        # print(r.text)
        url_viz = "http://cs207g10viz.us-east-1.elasticbeanstalk.com/visualization/demo/" + r.text + "/"
        print(url_viz)
        webbrowser.open_new(url_viz)
        return

# from chemkin_g10.chemkin import ReactionSystem
# if __name__ == '__main__':
#     T = 900
#     R = 8.314
#     concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
#     rsystem = ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
#     rsystem.buildFromXml("../tests/data/xml/rxns_reversible.xml", concs)
#     sim = Simulator(rsystem, 0.05)

#     sim.solveODE()
#     # sim.visualize_dev()
#     # sim.visualize()
#     # print(sim.yout)
#     # sim.plot_specie_all()
#     # sim.plot_reaction_all()
#     sim.plot_specie("H2O")
#     # print(sim.eq_diff)
#     # print(sim.check_equilibrium(5, 5e-11))
#     # print(sim.equilibrium_graph())
