# QuantumThermalMachines

This repository contains code for simulating and analyzing quantum thermal machines using the Julia programming language. Below you will find instructions on how to set up and run the code.

## Julia Installation

To run the code in this repository, you will need to have Julia installed on your system. If you haven't installed Julia yet, visit the [Julia website](https://julialang.org/) and follow the installation instructions provided.


Alternatively, you can use your package manager to install Julia. For example, on Debian based Linux Distros:
```
apt install julia
```

## Dependencies

Before running the code, you need to install the necessary dependencies. Open a terminal and launch the Julia REPL by typing `julia`. Once you're inside the REPL, activate the package manager by typing `]`. You should see a prompt like this:
```
(@v1.8) pkg> 
```

`@v1.8` indicates your version of Julia. To install the dependencies, enter the following command:
```
add Tables ITensors CSV DataFrames Roots LinearAlgebra
```

## Running the non-interactive system

From the repo directory type
```
julia code/NonInteractingQD.jl
```

The results will be saved in the results directory. To visualize the results, you can use the Jupyter Notebook provided the code directory. 


## Running the Interactive system

From the repo directory type
```
julia code/InteractingDSites.jl
```