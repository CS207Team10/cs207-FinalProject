# cs207 Team 10 Final Project: Chemical Kinetics library

[![Build Status](https://travis-ci.org/CS207Team10/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/CS207Team10/cs207-FinalProject)

[![Coverage Status](https://coveralls.io/repos/github/CS207Team10/cs207-FinalProject/badge.png?branch=master)](https://coveralls.io/github/CS207Team10/cs207-FinalProject?branch=master)

## Introduction:

This program is a chemical kinetics library, which can be used to calculate reaction rate coefficients, progress rate and reaction rate for a given system of chemical reactions. Our code is able to compute 3 different kinds of reaction rate coefficients: Constant reaction rate coefficients, Arrhenius reaction rate coefficients and Modified Arrhenius reaction rate coefficients as specified below. 

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align}&space;&k_{\textrm{const}}&space;=&space;k&space;\tag{constant}&space;\\&space;&k_{\textrm{arr}}&space;=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Arrhenius}&space;\\&space;&k_{\textrm{mod&space;arr}}&space;=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Modified&space;Arrhenius}&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align}&space;&k_{\textrm{const}}&space;=&space;k&space;\tag{constant}&space;\\&space;&k_{\textrm{arr}}&space;=&space;A&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Arrhenius}&space;\\&space;&k_{\textrm{mod&space;arr}}&space;=&space;A&space;T^{b}&space;\exp\left(-\frac{E}{RT}\right)&space;\tag{Modified&space;Arrhenius}&space;\end{align}" title="\begin{align} &k_{\textrm{const}} = k \tag{constant} \\ &k_{\textrm{arr}} = A \exp\left(-\frac{E}{RT}\right) \tag{Arrhenius} \\ &k_{\textrm{mod arr}} = A T^{b} \exp\left(-\frac{E}{RT}\right) \tag{Modified Arrhenius} \end{align}" /></a>

Each variable stands for A: Arrhenius prefactor, b: Modified Arrhenius parameter, E: Activation Energy, T: Temperature, and R: Ideal gas constant.

Progress rates and reaction rates of elementary reactions and irreversible reactions can also be handled by this program. 


A parse function is also included in our program, which was designed to parse an XML input file provided by users. The function handles the inputs by reading in all the needed data and store them for later calculation.


## Installation:

Our program includes two files: chemkin.py and test_chemkin.py. All the classes and functions necessary to run the library were stored in chemkin.py. All the tests were stored in test_chemkin.py. Both files can be downloaded from our repository. 

To run the program, first make sure you have Python installed on your machine. Then in your terminal, go to the directory where you store the two files. Run `python chemkin.py`. To run tests, run `pytest --cov test_chemkin.py`.


## Basic Usage and Examples: 


Our program includes three separate classes: `ChemUtil`, `Reaction` and `ReactionSystem`. 

### `ChemUtil`

`ChemUtil` is a class that contains necessary functions to compute all the coefficients and rates. It includes functions: `k_const`, `k_arr`, `k_mod_arr`, `progress_rate`, `reaction_rate` and `parse`.

### `Reaction`:

`Reaction` is a class that can be used to create a Reaction object. It includes functions: `updateCoeff` and `updateReaction`.

### `ReactionSystem`:

`ReactionSystem` is a class that represents a ReactionSystem object. It includes functions: `buildFromList`, `buildFromXml`, `getProgressRate` and `getReactionRate`.

### Example Output:

![alt text](https://github.com/CS207Team10/cs207-FinalProject/blob/master/example%20output.png)

