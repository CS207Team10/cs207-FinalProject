# Milestone 2

[Milestone 2 Rubric](https://github.com/IACS-CS-207/cs207-F17/blob/master/project/milestone2/milestone2_rubric.md)


## Proposed Feature
You should motivate and describe the feature
Explain how the feature will fit into your code base (and package)
Discuss the modules that you will write to realize your feature
Map out the methods you plan on implementing
Overview how you envision the user to use your new feature
Discuss any external dependencies that your feature will require


Based on Le Chatelier's principle, changing the concentration of a chemical will have an effect to the reaction equilibrium. Thus, the reaction rate, extent, and yield of products will be altered corresponding to the impact on the system. 

For irreversible reactions, the users may want to know the progress of the reaction by checking the concentration for each reactant and product. For reversible reactions, the users may also want to check the concentration for each specie, in order to determine whether the reaction has reached the equilibrium. 

First, our team plan to add a new feature to compute the concentration change for each specie at a specific time which will be provided by users. Since the concentration change of specie i is determined by the ODE: 

<a href="http://www.codecogs.com/eqnedit.php?latex=$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t}&space;=&space;f_{i}\left(\mathrm{x},&space;T\right),&space;\qquad&space;i&space;=&space;1,&space;\ldots,&space;N.$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t}&space;=&space;f_{i}\left(\mathrm{x},&space;T\right),&space;\qquad&space;i&space;=&space;1,&space;\ldots,&space;N.$$" title="$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t} = f_{i}\left(\mathrm{x}, T\right), \qquad i = 1, \ldots, N.$$" /></a>

and we already had the feature to compute the reaction rates, we can easily get the concentration change for each specie after a given time interval by intergrating the reaction rates by time.

Second, we also plan to add a new feature to determine whether a reversible reaction has reached its equilibrium. In chemistry, chemical equilibrium is the state in which both reactants and products are present in concentrations which have no further tendency to change with time. Usually, this state results when the forward reaction proceeds at the same rate as the reverse reaction. 

For example, a reversible reaction:

<a href="http://www.codecogs.com/eqnedit.php?latex=$$aA&space;&plus;&space;bB&space;\rightleftharpoons&space;cC&space;&plus;&space;dD$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$aA&space;&plus;&space;bB&space;\rightleftharpoons&space;cC&space;&plus;&space;dD$$" title="$$aA + bB \rightleftharpoons cC + dD$$" /></a>

<a href="http://www.codecogs.com/eqnedit.php?latex=$$\text{forward&space;progress&space;rate&space;=&space;}k^f&space;[A]^a[B]^b$$&space;\newline&space;$$\text{backward&space;progress&space;rate&space;=&space;}k^b&space;[C]^c[D]^d$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$\text{forward&space;progress&space;rate&space;=&space;}k^f&space;[A]^a[B]^b$$&space;\newline&space;$$&space\text{backward&space;progress&space;rate&space;=&space;}k^b&space;[C]^c[D]^d$$" title="$$\text{forward progress rate = }k^f [A]^a[B]^b$$ \newline $$\text{backward progress rate = }k^b [C]^c[D]^d$$" /></a>

When forward progress rate is equal to backward progress rate, the reaction reaches the equilibrium:

<a href="http://www.codecogs.com/eqnedit.php?latex=$$k^f&space;[A]^a[B]^b&space;=&space;k^b&space;[C]^c[D]^d$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$k^f&space;[A]^a[B]^b&space;=&space;k^b&space;[C]^c[D]^d$$" title="$$k^f [A]^a[B]^b = k^b [C]^c[D]^d$$" /></a>


Since our module already has the functions to compute forward reaction rate coefficient and backward reaction rate coefficient, if the users provide a time t, we can determine whether the reaction has reached equilibrium by first calculating the concentration for each specie at that time and then comparing the forward progress rate and the backward progress rate.

We will add a new class called ``equilibrium`` under  ``chemkin`` module. It will consist of several new methods, including one method to compute concentration change, which will take a time t as an input parameter, one method to determine equilibrium, and probably some other helper methods.

To use our new features, the users can first provide the reaction equation, initial concentration for each specie, necessary constants(A, b, E, T, R) and a time t, and then call the concentration change function to compute concentration change for each specie in this reaction, or directly call equilibrium function to see whether the reaction reaches the equilibrium at the given time t.



