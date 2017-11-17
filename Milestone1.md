# cs207 Team 10 Final Project: Chemical Kinetics library

[![Build Status](https://travis-ci.org/CS207Team10/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/CS207Team10/cs207-FinalProject)

[![Coverage Status](https://coveralls.io/repos/github/CS207Team10/cs207-FinalProject/badge.svg?branch=master&maxAge=1)](https://coveralls.io/github/CS207Team10/cs207-FinalProject?branch=master)

## Introduction:

This program is a chemical kinetics library, which can be used to calculate reaction rate coefficients, progress rate and reaction rate for a given system of chemical reactions. Our code is able to compute 3 different kinds of reaction rate coefficients: Constant reaction rate coefficients, Arrhenius reaction rate coefficients and Modified Arrhenius reaction rate coefficients as specified below. 

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align}&space;&k_{\textrm{const}}&space;=&space;k&space;\tag{constant}&space;\\&space;&k_{\textrm{arr}}&space;=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Arrhenius}&space;\\&space;&k_{\textrm{mod&space;arr}}&space;=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Modified&space;Arrhenius}&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;&k_{\textrm{const}}&space;=&space;k&space;\tag{constant}&space;\\&space;&k_{\textrm{arr}}&space;=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Arrhenius}&space;\\&space;&k_{\textrm{mod&space;arr}}&space;=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Modified&space;Arrhenius}&space;\end{align}" title="\begin{align} &k_{\textrm{const}} = k \tag{constant} \\ &k_{\textrm{arr}} = A \exp\left(-\frac{E}{RT}\right) \tag{Arrhenius} \\ &k_{\textrm{mod arr}} = A T^{b} \exp\left(-\frac{E}{RT}\right) \tag{Modified Arrhenius} \end{align}" /></a>

Each variable stands for A: Arrhenius prefactor, b: Modified Arrhenius parameter, E: Activation Energy, T: Temperature, and R: Ideal gas constant.

Progress rates and reaction rates of elementary reactions and irreversible reactions can also be handled by this program. 


A parse function is also included in our program, which was designed to parse an XML input file provided by users. The function handles the inputs by reading in all the needed data and store them for later calculation.


## Installation:

All the classes and functions necessary to run the library were stored in chemkin.py. All the tests were stored in test_chemkin.py. Both files can be downloaded from our repository, or you can obtain them here
[chemkin.py](https://raw.githubusercontent.com/CS207Team10/cs207-FinalProject/master/chemkin.py), [test_chemkin.py](https://raw.githubusercontent.com/CS207Team10/cs207-FinalProject/master/test_chemkin.py).

To use our program as library, move the two files into the same directory as your module, and
```python
import chemkin.py 
```
To run our test program, 
```
pytest --doctest-modules --cov --cov-report term-missing test_chemkin.py
```

## Basic Usage and Examples: 

Our program includes three separate classes: `ChemUtil`, `Reaction` and `ReactionSystem`. 

### ``ChemUtil`` class

`ChemUtil` is an util class that contains all necessary functions to compute the coefficients and rates. It includes functions: `k_const`, `k_arr`, `k_mod_arr`, `progress_rate`, `reaction_rate` and `parse`. 

### ``Reaction`` class

`Reaction` is a class that can be used to create a Reaction object. It includes functions: `updateCoeff` and `updateReaction`.

### ``ReactionSystem`` class

`ReactionSystem` is a class that represents a ReactionSystem object. It includes functions: `buildFromList`, `buildFromXml`, `getProgressRate` and `getReactionRate`.

### Examples

A typical workflow starts from initializing a `ReactionSystem` object. We need to set up all the needed variables: `T` (temperature), `R` (universal gas constant) and `concs` (the concentration of each species, and the order should be same with the one in the input file). Here's an example:

```python
import numpy as np
import chemkin as ck
T = 1500
R = 8.314
concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
rsystem = ck.ReactionSystem(T, R, concs)
```

Then, we feed an input `.xml` file (see [test1.xml](https://github.com/CS207Team10/cs207-FinalProject/blob/master/test1.xml) or [rxns.xml](https://github.com/CS207Team10/cs207-FinalProject/blob/master/rxns.xml) for the format) to the system by calling `buildFromXml` function:

```python
rsystem.buildFromXml("test1.xml")
```

We then can easily retrieve the progress rate and reaction rate by:

```python
rsystem.getProgressRate()
rsystem.getReactionRate()
```

If we want to get a summary of each reaction and also the whole system, we can do this by:

```python
print(rsystem)

```
```
The system:
2H2 + O2 =] 2OH + H2
	rate coeff: 70279405.1912
	rate coeff metadata: {'T': 1500, 'R': 8.314, 'type': 'modifiedArrhenius', 'A': 100000000.0, 'b': 0.5, 'E': 50000.0}
	reaction metadata: {'reversible': 'no', 'type': 'Elementary', 'id': 'reaction01'}
OH + HO2 =] H2O + O2
	rate coeff: 10000.0
	rate coeff metadata: {'T': 1500, 'R': 8.314, 'type': 'Constant'}
	reaction metadata: {'reversible': 'no', 'type': 'Elementary', 'id': 'reaction02'}
H2O + O2 =] HO2 + OH
	rate coeff: 4484938.47318
	rate coeff metadata: {'T': 1500, 'R': 8.314, 'type': 'Arrhenius', 'A': 10000000.0, 'E': 10000.0}
	reaction metadata: {'reversible': 'no', 'type': 'Elementary', 'id': 'reaction03'}
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
