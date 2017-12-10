# cs207 Team 10 Final Project: Chemical Kinetics library

[![Build Status](https://travis-ci.org/CS207Team10/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/CS207Team10/cs207-FinalProject)

[![Coverage Status](https://coveralls.io/repos/github/CS207Team10/cs207-FinalProject/badge.svg?branch=master&maxAge=1)](https://coveralls.io/github/CS207Team10/cs207-FinalProject?branch=master)


Introduction: Describe what problem the code is solving. You may borrow the Latex expressions from my lecture notes. Discuss in broad strokes what the purpose of the code is along with any features. Do not describe the details of the code yet.

Installation: Tell the user how to install the library. Be thorough here. If there is more than one way to find and use the code, then please clearly discuss each method.
Suggest a preferred installation method. Tell the user how they can contribute to the development version if they so desire. Tell the user how to run the test suite. Be explicit on any dependencies.

Basic Usage and Examples: Provide a few examples on using your software in some common situations. You may want to show how the code works with a small set of reactions. This is where you would introduce a toy demo problem and show off basic use cases. You would also highlight how to use your new feature here as well. If your code produces output or figures, you may want to include an examples/ directory in your repo so the user can try to reproduce the results in this section.

New Feature: This is where you discuss the final feature. This section should contain the following information:
Motivation of the new feature.
Implementation details including any new modules, classes, or methods. I can find the exact implementation in your code base, so your job here is to make sure I can understand why you made certain design decisions and how everything works together.


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

The library is able to compute 4 different kinds of reaction rate coefficient: Constant reaction rate coefficient, Arrhenius reaction rate coefficient, Modified Arrhenius reaction rate coefficient and backward reaction rate coefficient as specified below. 

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\begin{align*}&space;k_{\textrm{const}}&space;&=&space;k&space;&&&space;\text{constant}&space;\\&space;k_{\textrm{arr}}&space;&=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Arrhenius}&space;\\&space;k_{\textrm{mod&space;arr}}&space;&=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Modified&space;Arrhenius}&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;k_{\textrm{const}}&space;&=&space;k&space;&&&space;\text{constant}&space;\\&space;k_{\textrm{arr}}&space;&=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Arrhenius}&space;\\&space;k_{\textrm{mod&space;arr}}&space;&=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Modified&space;Arrhenius}&space;\end{align*}" title="\begin{align*} k_{\textrm{const}} &= k && \text{constant} \\ k_{\textrm{arr}} &= A \exp\left(-\frac{E}{RT}\right) && \text{Arrhenius} \\ k_{\textrm{mod arr}} &= A T^{b} \exp\left(-\frac{E}{RT}\right) && \text{Modified Arrhenius} \end{align*}" /></a>

Each variable stands for A: Arrhenius prefactor, b: Modified Arrhenius parameter, E: Activation Energy, T: Temperature, and R: Ideal gas constant.

For an elementary reaction (and only elementary reactions), we have 
<img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7B%5Cleft(b%5Cright)%7D%20%3D%20%5Cfrac%7Bk_%7Bj%7D%5E%7B%5Cleft(f%5Cright)%7D%7D%7Bk_%7Bj%7D%5E%7Be%7D%7D%2C%20%5Cqquad%20j%20%3D1%2C%20%5Cldots%2C%20M" alt="k_{j}^{\left(b\right)} = \frac{k_{j}^{\left(f\right)}}{k_{j}^{e}}, \qquad j =1, \ldots, M" />
    
where <img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7Be%7D" alt="k_{j}^{e}" /> is the *equilibrium coefficient* for reaction <img src="https://tex.s2cms.ru/svg/j" alt="j" />.
The final expression for the equilibrium coefficient is, 
   
<img src="https://tex.s2cms.ru/svg/k_%7Bj%7D%5E%7Be%7D%20%3D%20%5Cleft(%5Cfrac%7Bp_%7B0%7D%7D%7BRT%7D%5Cright)%5E%7B%5Cgamma_%7Bj%7D%7D%5Cexp%5Cleft(%5Cfrac%7B%5CDelta%20S_%7Bj%7D%7D%7BR%7D%20-%20%5Cfrac%7B%5CDelta%20H_%7Bj%7D%7D%7BRT%7D%5Cright)%2C%20%5Cqquad%20j%20%3D1%2C%20%5Cldots%2C%20M" alt="k_{j}^{e} = \left(\frac{p_{0}}{RT}\right)^{\gamma_{j}}\exp\left(\frac{\Delta S_{j}}{R} - \frac{\Delta H_{j}}{RT}\right), \qquad j =1, \ldots, M" />

where <img src="https://tex.s2cms.ru/svg/%5Cgamma_%7Bj%7D%20%3D%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%7B%5Cnu_%7Bij%7D%7D" alt="\gamma_{j} = \sum_{i=1}^{N}{\nu_{ij}}" /> and <img src="https://tex.s2cms.ru/svg/p_%7B0%7D" alt="p_{0}" /> is the pressure of the reactor (take it to be <img src="https://tex.s2cms.ru/svg/10%5E%7B5%7D" alt="10^{5}" /> Pa). <img src="https://tex.s2cms.ru/svg/%5C%5CDelta%20S_%7Bj%7D" alt="\\Delta S_{j}" /> is the entropy change of reaction <img src="https://tex.s2cms.ru/svg/j" alt="j" /> and <img src="https://tex.s2cms.ru/svg/%5C%5CDelta%20H_%7Bj%7D" alt="\\Delta H_{j}" /> the enthalpy change of reaction <img src="https://tex.s2cms.ru/svg/j" alt="j" />.


### Progress rate for irreversible and reversible reactions
 
 We used the principle of mass action to obtain the progress rate of each reaction.
 
In essence, we assert that the progress rate of a reaction is proportional to the concentrations of the reactants.
Thus, the forward progress rate is:

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

A parse function is also included in our library, which was designed to parse an XML input file provided by users. The function handles inputs by reading in all the needed data and store them for later calculation. 

Our program includes four separate modules: `chemkin`, `computation`, `db` and `thermo`. 

### ``chemkin`` class

`chemkin` is a class that represents a `ReactionSystem` object and a `Reaction` object. It includes functions: `updateCoeff`, `updateReaction`, `buildFromList`, `buildFromXml`, `getProgressRate`, `getReactionRate` and `parse`, which helps to create a `ReactionSystem` by parsing .xml input file.

### ``computation`` module

`computation` module contains all necessary functions to compute the coefficients and rates. It includes functions: `k_const`, `k_arr`, `k_mod_arr`, `progress_rate`, `reaction_rate`.

 
### ``thermo`` module

`thermo` includes functions: `H_over_RT`, `S_over_R`, `backward_coeffs`

### ```db``` module

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

   We can easily add a new method of computing reaction rate coefficient to our `ChemUtil` class.

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

### ODE Solver

### Plot

Based on Le Chatelier's principle, changing the concentration of a chemical will have an effect to the reaction equilibrium. Thus, the reaction rate, extent, and yield of products will be altered corresponding to the impact on the system. 

For irreversible reactions, users may want to know the progress of the reaction by checking the concentration for each reactant and product. For reversible reactions, users may also want to check the concentration for each specie, in order to determine whether the reaction has reached the equilibrium. 

First, our team plan to add a new feature to compute the concentration change for each specie after a time interval which will be provided by users. Since the concentration change of specie i is determined by the ODE: 

<a href="http://www.codecogs.com/eqnedit.php?latex=$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t}&space;=&space;f_{i}\left(\mathrm{x},&space;T\right),&space;\qquad&space;i&space;=&space;1,&space;\ldots,&space;N.$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t}&space;=&space;f_{i}\left(\mathrm{x},&space;T\right),&space;\qquad&space;i&space;=&space;1,&space;\ldots,&space;N.$$" title="$$\frac{\mathrm{d}x_{i}}{\mathrm{d}t} = f_{i}\left(\mathrm{x}, T\right), \qquad i = 1, \ldots, N.$$" /></a>

and we already had the feature to compute the reaction rates, we can easily get the concentration change for each specie after a given time interval by intergrating the reaction rates by time.

Second, we also plan to add a new feature to determine whether a reversible reaction has reached its equilibrium. In chemistry, chemical equilibrium is the state in which both reactants and products are present in concentrations which have no further tendency to change with time. Usually, this state results when the forward reaction proceeds at the same rate as the reverse reaction. 

For example, a reversible reaction:

<a href="http://www.codecogs.com/eqnedit.php?latex=$$aA&space;&plus;&space;bB&space;\rightleftharpoons&space;cC&space;&plus;&space;dD$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$aA&space;&plus;&space;bB&space;\rightleftharpoons&space;cC&space;&plus;&space;dD$$" title="$$aA + bB \rightleftharpoons cC + dD$$" /></a>

<a href="http://www.codecogs.com/eqnedit.php?latex=$$\text{forward&space;progress&space;rate&space;=&space;}k^f&space;[A]^a[B]^b$$&space;\newline&space;$$\text{backward&space;progress&space;rate&space;=&space;}k^b&space;[C]^c[D]^d$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$\text{forward&space;progress&space;rate&space;=&space;}k^f&space;[A]^a[B]^b$$&space;\newline&space;$$&space\text{backward&space;progress&space;rate&space;=&space;}k^b&space;[C]^c[D]^d$$" title="$$\text{forward progress rate = }k^f [A]^a[B]^b$$ \newline $$\text{backward progress rate = }k^b [C]^c[D]^d$$" /></a>

When forward progress rate is equal to backward progress rate, the reaction reaches the equilibrium:

<a href="http://www.codecogs.com/eqnedit.php?latex=$$k^f&space;[A]^a[B]^b&space;=&space;k^b&space;[C]^c[D]^d$$" target="_blank"><img src="http://latex.codecogs.com/gif.latex?$$k^f&space;[A]^a[B]^b&space;=&space;k^b&space;[C]^c[D]^d$$" title="$$k^f [A]^a[B]^b = k^b [C]^c[D]^d$$" /></a>


Since our module already has the functions to compute forward reaction rate coefficient and backward reaction rate coefficient, if the users provide a time t, we can determine whether the reaction has reached equilibrium by first calculating the concentration for each specie at that time and then comparing the forward progress rate and the backward progress rate.

We will add a new module called ``equilibrium``. It will consist of several new methods, including one method to compute concentration change, which will take a time t as an input parameter, one method to determine equilibrium, and probably some other helper methods.

To use our new features, the users can first provide the reaction equation, initial concentration for each specie, necessary constants(A, b, E, T, R) and a time t, and then call the concentration change function to compute concentration change for each specie in this reaction, or directly call equilibrium function to see whether the reaction reaches the equilibrium at the given time t.

### equilibrium
It is useful to know if the system has reached equilibrium at the end of the simulation at time "t".
The new feature named "equilibrium" judges if the sytem has reached equilibrium.
Here, we consider a simulation with time-interval (0,t).
If the largest concentration among chemical species at time "t" is "C", then the characteristic slope of the c(t) curves can be calculated as "C/t".
We judge the system to be in equilibrium if all the slopes of the concentrations at the last two time steps are less than the critical slope value "1e-8xC/t".

### Web?
