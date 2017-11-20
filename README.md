# cs207 Team 10 Final Project: Chemical Kinetics library

[![Build Status](https://travis-ci.org/CS207Team10/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/CS207Team10/cs207-FinalProject)

[![Coverage Status](https://coveralls.io/repos/github/CS207Team10/cs207-FinalProject/badge.svg?branch=master&maxAge=1)](https://coveralls.io/github/CS207Team10/cs207-FinalProject?branch=master)


## Introduction:

This program is a chemical kinetics library, which can be used to calculate reaction rate coefficients, progress rate and reaction rate for a given system of chemical reactions. 

The functionalities include computing 3 different kinds of reaction rate coefficients: Constant reaction rate coefficients, Arrhenius reaction rate coefficients and Modified Arrhenius reaction rate coefficients as specified below. 

<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\begin{align*}&space;k_{\textrm{const}}&space;&=&space;k&space;&&&space;\text{constant}&space;\\&space;k_{\textrm{arr}}&space;&=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Arrhenius}&space;\\&space;k_{\textrm{mod&space;arr}}&space;&=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Modified&space;Arrhenius}&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;\begin{align*}&space;k_{\textrm{const}}&space;&=&space;k&space;&&&space;\text{constant}&space;\\&space;k_{\textrm{arr}}&space;&=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Arrhenius}&space;\\&space;k_{\textrm{mod&space;arr}}&space;&=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;&&&space;\text{Modified&space;Arrhenius}&space;\end{align*}" title="\begin{align*} k_{\textrm{const}} &= k && \text{constant} \\ k_{\textrm{arr}} &= A \exp\left(-\frac{E}{RT}\right) && \text{Arrhenius} \\ k_{\textrm{mod arr}} &= A T^{b} \exp\left(-\frac{E}{RT}\right) && \text{Modified Arrhenius} \end{align*}" /></a>

Each variable stands for A: Arrhenius prefactor, b: Modified Arrhenius parameter, E: Activation Energy, T: Temperature, and R: Ideal gas constant.

Progress rates and reaction rates of both **irreversible** and **reversible** reactions can also be computed by the library. Specifically, for **reversible** function, the total progress rate is coumputed by
<img src="https://github.com/CS207Team10/cs207-FinalProject/blob/master/docs/images/backward.svg">

And we're using the NASA Polynomial Coefficients for computing and entropy and enthalpy change of each reaction.

A parse function is also included in our library, which was designed to parse an XML input file provided by users. The function handles inputs by reading in all the needed data and store them for later calculation. 


## Installation:

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

## Basic Usage and Examples: 

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

## Proposed Feature

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


