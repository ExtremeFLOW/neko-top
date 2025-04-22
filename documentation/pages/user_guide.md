# User guide {#user_guide}

Neko-TOP is an extension of the Neko library, and as such it requires the Neko
library to be installed. The Neko library also has a number
of dependencies, such as the JSON-Fortran library, which is used for reading and
writing JSON files, and the PFUnit library, which is used for unit testing.

- \subpage installation
  - [Dependencies](@ref installation-external)
  - [Quick-start compilation](@ref installation-quick)
- \subpage examples
  - [Running the examples](@ref examples-running)
  - [Creating a new case](@ref examples-new)
  - [Case files and meshes](@ref examples-data)
  - [Advanced example setups](@ref examples-advanced)
- \subpage clusters

\note the remaining structure here should be discussed. I suggest we follow
in suit with `neko`, where they essentially break up the user guide based on
the top level structure of the case file. In any case, this would be where we
introduce these subpages
\subpage mapping_cascade
\subpage objectives_and_constraints
\subpage simulation_components
