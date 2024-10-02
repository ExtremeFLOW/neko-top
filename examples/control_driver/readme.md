# Proof of concept for a steady state optimization problem

Here we follow the classic "rugby ball" example by Borrvall & Petersson
https://doi.org/10.1002/fld.426



This is example is primarily used to summarize the list of remaining tasks.

However, the source code is littered with TODO comments and further details.

All these task have been labeled from 
	one star (*), implying an easy fix

	three stars(***), implying something more challenging

## Todo (in Neko):


 Simcomp_handler (like Tim made for source terms)

Potentially steady states


## Todo (in Neko-top):

Masks

incorporate the PDE filters


incorporate different mappings

Make the "problem" module

Handle mapping to multiple coefficients 

Handle multiple constraints (such that we have multiple sensitivities)

log file write with KKT, constraint values etc

merge in the passive scalars
