# PyBofaa

## Description

An implementatoin of the multi-objective metaheuristcs ***Multi-objective Genetic Algorithm (NSGA-II)*** and ***Multi-objective Bat Algorithm (MOBA)***, to the *DNA fragment assembly problem - FAP*, while taking three objectives: Minimizing the number of oweverlaping distant fragments, Maximizing the number of oweverlaping adjacent fragments, and Minimizing the number of contigs; in order to solve the FAP problem and optimize the solution.

## Tools

* Language: Python 3.6
* Libraries: typing, math, random
* Compiler: PyPy, version 7.3.2

For perofrmance reasons we used the PyPy compiler, that reduced the execution time of the code during the tests on the listed benchmarks, to 7x times less, comparing to the execution time when using the standard Python3 compiler.

## How To Read The Code

* `config.py`
  > Contains all the parameters values for both the ***NSGA-II*** and ***MOBA*** algorithm, the list of all the benchmark files, and the selected file for the benchmark.

* `run_nsga2.py`
  > For compiling the ***NSGA-II*** algorithm

* `run_bat_algorithm.py`
  > For compiling the ***MOBA*** algorithm

* `models/`
  * `Fragment.py`
    > A class for modeling the DNA fragment that is stored as a *String*, in which all the possible related data is calculated and stored, to make it for debbuging and code tracking.

  * `Solution.py`
    > A class for modeling a possible solution to our problem, that is inintialy a possible sequencing of our DNA fragments, the solution is a ***list of integers***, where each element represents the index of the corresponging fragment. All the relative data is stored into the object.

* `use/`
  > Contains the `scoring.py` and `tools.py`, contains all the neccesary algorithms ro read data, and calculate the overlap scores.
  
* `algorithm/`
  * `MultiObjective.py`
    > A class with all the methods needed for the multi objectivity problems, used by both the ***NSGA-II*** and ***MOBA*** algorithms.

  * `BatAlgorithm.py`, `NsGa2.py`
    > Two classes, that each one contains the required methods for implementing the two algorithms.

* `benchmarks`
  > The file benchmarks, used to test the algorithm.

## Run Tests

To run tests, just use the two files `run_nsga2.py` for using the ***NSGA-II*** algorithm, and  `run_bat_algorithm.py` for using the ***MOBA*** algorithm.

You can simply change the `confi.py` file format, but keep the execution order of the functions, since the `run_nsga2.py`, `run_bat_algorithm.py` represents the algorithms.

## Contributors

[Herhar Fares](https://github.com/HerharFares)  
[Brahim Benyoub](mailto:brahim.weldtest@gmail.com)

## Note

***Use the PyPy compiler, because it will take a very long time for the code to be executed, using the standard Python compiler.***
