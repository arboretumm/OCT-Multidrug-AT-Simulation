# Code

This is a file of matlab code that runs ODE simulations of different two-drug regimens.

## Installation

Download all the folders and add to a directory that matlab can access.

## Directory
base run function - contains all the functions that can run a single regimen simulation. Needed for running multiple simulation runs. 

equation functions - contains the ODE system equations used for different drugging approaches

event functions - contains the functions that define different event triggers for the code. eventSubpopulationLimit define an event based off of the size of a subpopulation and the folders in concentration triggers define an event based off of the concentration of a drug.

other figure gen - contains the code to make a bar chart of the different population initial conditions.

other tools - contains the code for practical versions of regimens, for the standard-of-care regimens (referred to as maximalTreatment here), for the theoretical best and worst, and for comparing the outcomes of different regimen simulations (OCT shape comparison folder)

population parameters - contains the code for generation a population struct

simulation runs - contains various multi and mega runs of the simulation. multiruns are defined as runs containing multiple runs of a single base function over different parameters, while megaruns contain multiple multiruns. There is some unhygienic naming here - both multirun and multisimrun were used to refer to multiruns. 

## License

[MIT](https://choosealicense.com/licenses/mit/)
