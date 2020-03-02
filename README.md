# parallel_tempering

C implementation of parallel tempering for the 2D Ferromagnetic Ising Model.  Acceptance probabilities, A(T), functions of temperature, T, are estimated given a temperature set.  ```results``` includes example ```config``` and ```output``` files for geometric and optimized temperature sets, generating a results similar to Fig. 5 in [Katzgraber, et al., 2006](https://arxiv.org/abs/cond-mat/0602085): 
![fig5_production.png](results/fig5_production.png) 

## Requirements

- gcc 8.1.0
- GSL
- openMP

## Installation

Clone the repo to a local directory.

## Usage

Edit system size, temperature set, and other simulation parameters in ```config.txt```.  In a bash shell, run ```make``` and then ```./parallel_tempering```. 

