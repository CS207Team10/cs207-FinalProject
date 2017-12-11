[![Build Status](https://travis-ci.org/CS207Team10/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/CS207Team10/cs207-FinalProject)

[![Coverage Status](https://coveralls.io/repos/github/CS207Team10/cs207-FinalProject/badge.svg?branch=master&maxAge=1)](https://coveralls.io/github/CS207Team10/cs207-FinalProject?branch=master)

# cs207 Team 10 Final Project: Chemical Kinetics library

#### Team Members:  Hidenori Tanaka, Jiachen Song and Xiangru Shu

## Introduction

This program is a chemical kinetics library, which could be easily installed by users and used for various applications.

Basic usage includes:

 * Constant reaction rate coefficient
 * Arrhenius reaction rate coefficient
 * Modified Arrhenius reaction rate coefficient
 * Backward reaction rate coefficient
 * Progress rate for irreversible and reversible reactions
 * Reaction rate for species
 
 ### Reaction rate coefficients

The library is able to compute 4 different kinds of reaction rate coefficient: constant reaction rate coefficient, arrhenius reaction rate coefficient, modified Arrhenius reaction rate coefficient and backward reaction rate coefficient as specified below. 

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\begin{align*}&space;k_{\textrm{const}}&space;&=&space;k&space;&&&space;\text{constant}&space;\\&space;k_{\textrm{arr}}&space;&=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Arrhenius}&space;\\&space;k_{\textrm{mod&space;arr}}&space;&=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Modified&space;Arrhenius}&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;k_{\textrm{const}}&space;&=&space;k&space;&&&space;\text{constant}&space;\\&space;k_{\textrm{arr}}&space;&=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Arrhenius}&space;\\&space;k_{\textrm{mod&space;arr}}&space;&=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Modified&space;Arrhenius}&space;\end{align*}" title="\begin{align*} k_{\textrm{const}} &= k && \text{constant} \\ k_{\textrm{arr}} &= A \exp\left(-\frac{E}{RT}\right) && \text{Arrhenius} \\ k_{\textrm{mod arr}} &= A T^{b} \exp\left(-\frac{E}{RT}\right) && \text{Modified Arrhenius} \end{align*}" /></a>

Each variable stands for A: Arrhenius prefactor, b: Modified Arrhenius parameter, E: Activation Energy, T: Temperature, and R: Ideal gas constant.

For an elementary reaction, we have 
<img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7B%5Cleft(b%5Cright)%7D%20%3D%20%5Cfrac%7Bk_%7Bj%7D%5E%7B%5Cleft(f%5Cright)%7D%7D%7Bk_%7Bj%7D%5E%7Be%7D%7D%2C%20%5Cqquad%20j%20%3D1%2C%20%5Cldots%2C%20M" alt="k_{j}^{\left(b\right)} = \frac{k_{j}^{\left(f\right)}}{k_{j}^{e}}, \qquad j =1, \ldots, M" />
    
where <img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7Be%7D" alt="k_{j}^{e}" /> is the *equilibrium coefficient* for reaction <img src="https://tex.s2cms.ru/svg/j" alt="j" />.
The final expression for the equilibrium coefficient is, 
   
<img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7Be%7D%20%3D%20%5Cleft(%5Cfrac%7Bp_%7B0%7D%7D%7BRT%7D%5Cright)%5E%7B%5Cgamma_%7Bj%7D%7D%5Cexp%5Cleft(%5Cfrac%7B%5CDelta%20S_%7Bj%7D%7D%7BR%7D%20-%20%5Cfrac%7B%5CDelta%20H_%7Bj%7D%7D%7BRT%7D%5Cright)%2C%20%5Cqquad%20j%20%3D1%2C%20%5Cldots%2C%20M" alt="k_{j}^{e} = \left(\frac{p_{0}}{RT}\right)^{\gamma_{j}}\exp\left(\frac{\Delta S_{j}}{R} - \frac{\Delta H_{j}}{RT}\right), \qquad j =1, \ldots, M" />

where <img src="https://tex.s2cms.ru/svg/%5Cgamma_%7Bj%7D%20%3D%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%7B%5Cnu_%7Bij%7D%7D" alt="\gamma_{j} = \sum_{i=1}^{N}{\nu_{ij}}" /> and <img src="https://tex.s2cms.ru/svg/p_%7B0%7D" alt="p_{0}" /> is the pressure of the reactor (take it to be <img src="https://tex.s2cms.ru/svg/10%5E%7B5%7D" alt="10^{5}" /> Pa). <img src="https://tex.s2cms.ru/svg/%5CDelta%20S_%7Bj%7D" alt="\Delta S_{j}" /> is the entropy change of reaction <img src="https://tex.s2cms.ru/svg/j" alt="j" /> and <img src="https://tex.s2cms.ru/svg/%5CDelta%20H_%7Bj%7D" alt="\Delta H_{j}" /> the enthalpy change of reaction <img src="https://tex.s2cms.ru/svg/j" alt="j" />.


### Progress rate for irreversible and reversible reactions
 
 We used the principle of mass action to obtain the progress rate of each reaction.
 
Basically, the progress rate of a reaction is proportional to the concentrations of the reactants. Thus, the forward progress rate is:

 <img src="https://tex.s2cms.ru/svg/r_%7Bj%7D%20%3D%20k_%7Bj%7D%5E%7B(f)%7D%5Cprod_%7Bi%3D1%7D%5E%7BN%7D%7Bx_%7Bi%7D%5E%7B%5Cnu_%7Bij%7D%5E%7B%5Cprime%7D%7D%7D%2C%20%5Cqquad%20j%20%3D%201%2C%5Cldots%2C%20M" alt="r_{j} = k_{j}^{(f)}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime}}}, \qquad j = 1,\ldots, M" />

where <img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7B(f)%7D" alt="k_{j}^{(f)}" /> is the forward reaction rate coefficient.

For irreverisible reactions, the total progress rate is equal to the forward progress rate. But in reality, it is often the case that the products can react and produce the reactants. This is called a reversible reaction.

For reversible reaction, the total progress rate is:

<img src="https://tex.s2cms.ru/svg/r_%7Bj%7D%20%3D%20k_%7Bj%7D%5E%7B%5Cleft(f%5Cright)%7D%5Cprod_%7Bi%3D1%7D%5E%7BN%7D%7Bx_%7Bi%7D%5E%7B%5Cnu_%7Bij%7D%5E%7B%5Cprime%7D%7D%7D%20-%20k_%7Bj%7D%5E%7B%5Cleft(b%5Cright)%7D%5Cprod_%7Bi%3D1%7D%5E%7BN%7D%7Bx_%7Bi%7D%5E%7B%5Cnu_%7Bij%7D%5E%7B%5Cprime%5Cprime%7D%7D%7D%2C%20%5Cqquad%20j%20%3D%201%2C%5Cldots%2C%20M." alt="r_{j} = k_{j}^{\left(f\right)}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime}}} - k_{j}^{\left(b\right)}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime\prime}}}, \qquad j = 1,\ldots, M." />

where <img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7B%5Cleft(b%5Cright)%7D" alt="k_{j}^{\left(b\right)}" /> is the backward reaction rate coefficient.

 ### Reaction rate for species
 
The reaction rate for each species was given as a linear combination of reaction progress rates.  

That is, 

<img src="https://tex.s2cms.ru/svg/f_%7Bi%7D%5Cleft(%5Cmathbf%7Bx%7D%2C%20T%5Cright)%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7BM%7D%7B%5Cnu_%7Bij%7Dr_%7Bj%7D%7D%2C%20%5Cqquad%20i%20%3D%201%5Cldots%2C%20N" alt="f_{i}\left(\mathbf{x}, T\right) = \sum_{j=1}^{M}{\nu_{ij}r_{j}}, \qquad i = 1\ldots, N" /> 
   
where <img src="https://tex.s2cms.ru/svg/%5Cnu_%7Bij%7D%20%3D%20%5Cnu_%7Bij%7D%5E%7B%5Cprime%5Cprime%7D%20-%20%5Cnu_%7Bij%7D%5E%7B%5Cprime%7D" alt="\nu_{ij} = \nu_{ij}^{\prime\prime} - \nu_{ij}^{\prime}" />.



## Installation

All the classes and functions necessary to run the library were stored in chemkin_g10. All the tests were stored in tests folder. 

To install program as library (make sure you have pip installed for python3.5+)
```
pip install chemkin_g10 
```
Two samples were provided in the **samples** directory: irreversible.py and reversible.py.

To run our samples, go to the **samples** directory, and run
```
python irreversible.py
```
```
python reversible.py
```
To run our test program, simply run
```
pytest
```

## Basic Usage and Examples

Our library includes four separate modules: `chemkin`, `computation`, `db` and `thermo`. 

### chemkin module

`chemkin` is a module that contains a `ReactionSystem` class, a `Reaction` class and a `Simulator` class. 

1. `ReactionSystem` class:

    * `buildFromList`
    * `buildFromXml`
    * `getProgressRate`
    * `getReactionRate`
    * `parse`
    
2. `Reaction` class:

    * `updateCoeff`
    * `updateReaction`
    
3. `Simulator` class:
    * `solveODE`
    * `check_equilibrium`
    * `equilibrium_graph`
    * `plot_specie`
    * `plot_specie_all`
    * `plot_reaction_all`

### computation module

`computation` module contains all necessary functions to compute the coefficients and rates. 

It includes functions: 

* `k_const`
* `k_arr`
* `k_mod_arr`
* `progress_rate`
* `reaction_rate`
* `equilibrium_constant`

### thermo module

`thermo` includes functions: `H_over_RT`, `S_over_R`, `backward_coeffs`

### db module

```db``` contains the functions to read NASA polynomials from the database.

### Examples

A typical workflow starts from initializing a `ReactionSystem` object. We need to set up all the needed variables: `T` (temperature), `R` (universal gas constant) and `concs` (the concentration of each species, and the order should be same with the one in the input file). Here's an example:

```python
from chemkin_g10 import chemkin as ck
import numpy as np

T = 1500
R = 8.314
concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
rsystem = ck.ReactionSystem(T, R, "tests/data/db/nasa.sqlite") # your path to db file
```

Then, we feed an input `.xml` file (see [example1](https://github.com/CS207Team10/cs207-FinalProject/blob/master/tests/data/xml/rxns_short_1.xml) or [example2](https://github.com/CS207Team10/cs207-FinalProject/blob/master/tests/data/xml/rxns_reversible.xml) for the format) **and** a numpy array of concerntrations of each species to the system by calling `buildFromXml` function:

```python
rsystem.buildFromXml("tests/data/xml/rxns_reversible.xml", concs) # your path to xml file
```

We then can easily retrieve the progress rate and reaction rate by:

```python
print("Progress rate: \n", rsystem.getProgressRate(), "\n")
print("Reaction rate: \n", rsystem.getReactionRate(), "\n")
```

If we want to get a summary of each reaction and also the whole system, we can do this by:

```python
print("System info: \b", rsystem, "\n")
```

### Extensibility

1. Reversible/Non-elementary reactions
   
   We store these types of information in the metadata of each reaction, for example:  
   ```python
   rsystem.reactionList[0].reactMeta
   ```
   ```
   {'reversible': 'no', 'type': 'Elementary', 'id': 'reaction01'}
   ```

2. Reaction rate coefficients not discussed in class

   We can easily add a new method of computing reaction rate coefficient to our `computation` module.

3. Other extensions

   We can update the property (metadata) of a single reaction or its reaction rate coefficient (re-compute the coefficient), with any valid parameters. Here's an example of changing the coefficient of the 3rd reaction to be the same as the 1st reaction:
   
   ```python
   print(rsystem.reactionList[2].k, rsystem.reactionList[0].k)
   # using the parameters of the first equation
   rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
   # both ks are same now
   print(rsystem.reactionList[2].k, rsystem.reactionList[0].k) 
   ```
   ```
   4484938.47318 70279405.1912
   70279405.1912 70279405.1912
   ```

   Then rebuild the system with the new variables: 
   ```python
   rsystem.buildFromList(rsystem.reactionList)
   ```

## New Feature

We implemented a new class `Simulator` under the `chemkin` module. The class includes an ODE solver function, two equilibrium check functions and three plot functions.

To initialize a `Simulator` object, users should pass in a `ReactionSystem` object, a time stamp, a number of integration time steps(default = 100), a time scale(default = 1e9) and an equilibrium threshold(default = 1e-05). 

For example:
```
sim = Simulator(rsystem, 0.1, numSample=100, timeScale=1e9, eqThreshold=1e-05)
```

### ODE Solver

Based on Le Chatelier's principle, changing the concentration of a chemical will have an effect to the reaction equilibrium. Thus, the reaction rate, extent, and yield of products will be altered corresponding to the impact on the system. 

For irreversible reactions, users may want to know the progress of the reaction by checking the concentration for each reactant and product. For reversible reactions, users may also want to check the concentration for each specie, in order to determine whether the reaction has reached the equilibrium. 

Thus, we added a new feature to compute the concentration for each specie at each time stamp between time interval (0,t), and t will be entered by users. 

To implement this feature, we added a new function `solveODE` under the `Simulator` class. It has a nested function `fun`, which is to compute reaction rates. Then we used `odeint` function from library `scipy` to integrate reaction rates and then get concentrations for each specie at each time stamp. 

The formula of reaction rate for each specie is specified below:

<a href="http://www.codecogs.com/eqnedit.php?latex=$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t}&space;=&space;f_{i}\left(\mathrm{x},&space;T\right),&space;\qquad&space;i&space;=&space;1,&space;\ldots,&space;N.$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t}&space;=&space;f_{i}\left(\mathrm{x},&space;T\right),&space;\qquad&space;i&space;=&space;1,&space;\ldots,&space;N.$$" title="$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t} = f_{i}\left(\mathrm{x}, T\right), \qquad i = 1, \ldots, N.$$" /></a>

To use this function, users can simply call:

```
sim.solveODE()
```

Users can also get all the concentrations at each time stamp by calling the following methods:

```
print(sim.yout)
```

```
...
[  1.62865785e+00   9.55229994e-02   9.03020740e-03   2.06975115e-01
    1.22418085e+00   3.35632969e-01   4.49199685e-15   1.58667651e-20]
 [  1.62865785e+00   9.55229994e-02   9.03020740e-03   2.06975115e-01
    1.22418085e+00   3.35632969e-01   4.49199685e-15   1.58667651e-20]
 [  1.62865785e+00   9.55229994e-02   9.03020740e-03   2.06975115e-01
    1.22418085e+00   3.35632969e-01   4.49199685e-15   1.58667651e-20]]
```

The last array of concentrations is the final concentrations for each specie at time t.

### Equilibrium

In chemistry, chemical equilibrium is the state in which both reactants and products are present in concentrations which have no further tendency to change with time. Usually, this state results when the forward reaction proceeds at the same rate as the reverse reaction. 

To check equilibrium, we developed two methods. The first method is based on the definition of equilibrium and compare reaction quotient with equilibrium constants.

For example, a reversible reaction:

<a href="http://www.codecogs.com/eqnedit.php?latex=$$aA&space;&plus;&space;bB&space;\rightleftharpoons&space;cC&space;&plus;&space;dD$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$aA&space;&plus;&space;bB&space;\rightleftharpoons&space;cC&space;&plus;&space;dD$$" title="$$aA + bB \rightleftharpoons cC + dD$$" /></a>

The equation for reaction quotient is written by multiplying the concentrations for the species of the products and dividing by the concentrations of the reactants. If any component in the reaction has a coefficient, indicated above with lower case letters, the concentration is raised to the power of the coefficient. The above equation is therefore: 

<img src="https://tex.s2cms.ru/svg/Qc%3D%5Cfrac%7B%5BC%5D%5E%7Bc%7D%5BD%5D%5E%7Bd%7D%7D%7B%5BA%5D%5Ea%5BB%5D%5Eb%7D" alt="Qc=\frac{[C]^{c}[D]^{d}}{[A]^a[B]^b}" />

A comparison of reaction quotient(Q) with equilibrium constant(K) indicates which way the reaction shifts and which side of the reaction is favored:

* If  Q>K, then the reaction favors the reactants. 
* If  Q<K, then the reaction favors the products. 
* If  Q=K, then the reaction is already at equilibrium. There is no tendency to form more reactants or more products at this point. No side is favored and no shift occurs.

First, we added a new function `equilibrium_constant` under the `computation.py` module to compute the equilibrium constant(k) for each reaction.

Then, in the `solveODE` function, after calling `odeint` to integrate concentrations, we computed the reaction quotients for each reaction, and compare the reaction quotients with the equilibrium constants by computing the following percentage:

<img src="https://tex.s2cms.ru/svg/%5Cfrac%7B%5Ctext%7Breaction%20quotient%7D%20-%20%5Ctext%7Bequilibrium%20constant%7D%7D%7B%5Ctext%7Bequilibrium%20constant%7D%7D" alt="\frac{\text{reaction quotient} - \text{equilibrium constant}}{\text{equilibrium constant}}" />

If the percentage is less than the equilibrium threshold, then we consider the reaction has reached equilibrium. 

For each reaction, we found the first time stamp that the reaction reaches equilibrium and stored them into an array for future use.

Users can simply call the following methods to see the equilibrium time stamp for each reaction:
```
print(sim.eq_point)
```

```
[3.4343434343434346e-11, 6.7676767676767674e-11, 6.8686868686868691e-11, 6.2626262626262625e-11, 6.6666666666666669e-11, 5.858585858585858e-11, 5.858585858585858e-11, 6.4646464646464647e-11, 5.4545454545454549e-11, 6.5656565656565652e-11, 5.6565656565656564e-11]
```

Users may also want to check equilibrium with arbitrary time stamp, so we added a new function `check_equilibrium` which takes an index and a time stamp as parameters. By specifying index and time stamp, users can check equilibrium for a specific reaction at a specific time stamp. The function will simply compare this time stamp with the first equilibrium time stamp for that reaction and determine whether the reaction has reached equilibrium.

An example of checking equilibrium for reaction 5 at time 5e-11(5e-11 is less than 5.858585858585858e-11, so it has not reached equilibrium):
```
print(sim.check_equilibrium(5, 5e-11))
```

```
False
```

Furthermore, we came up with another method to check equilibrium, which was developed based on slopes of concentration plots.

If the largest concentration among chemical species at time t is C, then the characteristic slope of the c(t) curves can be calculated as C/t. We judge the system to be in equilibrium if all the slopes of the concentrations at the last two time steps are less than the critical slope value "1e-7*C/t". The choice of "1e-7" is our definition of equilibrium and the number could be changed to another small number.

Example:
```
print(sim.equilibrium_graph())
```
```
True
```

### Plot

It is always helpful for users to see how concentrations change over time graphically. Therefore, we added plot functions to visualize concentration change and we used `matplotlib` library to plot graphs.

To plot concentrations for the entire reaction system, we added function `plot_specie_all`, which will plot the concentrations for all the species in the system over time.
```
sim.plot_specie_all()
```

![plot_specie_all](https://github.com/CS207Team10/cs207-FinalProject/blob/master/images/Figure_1.png)

To plot concentration for an individual specie, we added function `plot_specie`, which will take an integer as parameter to specify which specie to plot. 

```
sim.plot_specie(4)
```

![plot_specie](https://github.com/CS207Team10/cs207-FinalProject/blob/master/images/Figure_2.png)

We also added a function `plot_reaction_all` to plot (reaction quotients - equilibrium constants)/ equilibrium constants.

```
sim.plot_reaction_all()
```

![plot_reaction_all](https://github.com/CS207Team10/cs207-FinalProject/blob/master/images/Figure3.png)

### Web?
