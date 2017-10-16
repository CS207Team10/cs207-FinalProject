# cs207-FinalProject: Chemical Kinetics library

[![Build Status](https://travis-ci.org/CS207Team10/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/CS207Team10/cs207-FinalProject)

[![Coverage Status](https://coveralls.io/repos/github/CS207Team10/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/CS207Team10/cs207-FinalProject?branch=master)

Rubric: [here](https://github.com/IACS-CS-207/cs207-F17/blob/master/project/milestone1_rubric.md)

## Introduction:
(Describe what problem the code is solving. You may borrow the Latex expressions from my lecture notes. Discuss in broad strokes what the purpose of the code is along with any features. Do not describe the details of the code yet.)

This program is a chemical kinetics library, which can be used to calculate reaction rate coefficients, progress rate and reaction rate for a given system of chemical reactions. Our code is able to compute 3 different kinds of reaction rate coefficients: Constant reaction rate coefficients, Arrhenius reaction rate coefficients and Modified Arrhenius reaction rate coefficients. Progress rates and reaction rates of elementary reactions and irreversible reactions can also be handle by this program. 

A parse function is also included in our program, which was designed to parse an XML input file provided by users. The function handles the inputs by reading in all the needed data and store them for later calculation.


## Installation:
(Describe where the code can be found and downloaded. Tell the user how to run the test suite. We are not releasing this code as a package yet, but when we do that this section will include instructions how how to install the package.)

Our program includes two files: chemkin.py and test_chemkin.py. All the classes and functions necessary to run the library were stored in chemkin.py. All the tests were stored in test_chemkin.py. Both files can be downloaded from our repository. 

To run the program, first make sure you have Python installed on your machine. Then in your terminal, go to the directory where you store the two files. Run `python chemkin.py`. To run tests, run `pytest --cov test_chemkin.py`.


## Basic Usage and Examples: 
(Provide a few examples on using your software in some common situations. You may want to show how the code works with a small set of reactions.)

Our program includes three separate classes: `ChemUtil`, `Reaction` and `ReactionSystem`. 

### `ChemUtil`

`ChemUtil` is a class that contains necessary functions to compute all the coefficients and rates. It includes functions: `k_const`, `k_arr`, `k_mod_arr`, `progress_rate`, `reaction_rate` and `parse`.

### `Reaction`:

`Reaction` is a class that can be used to create a Reaction object. It includes functions: `updateCoeff` and `updateReaction`.

### `ReactionSystem`:

`ReactionSystem` is a class that represents a ReactionSystem object. It includes functions: `buildFromList`, `buildFromXml`, `getProgressRate` and `getReactionRate`.




