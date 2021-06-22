# Dispatch Optimization with Pyomo

We cast the dispatch optimization as a mixed integer linear programming problem and solve it using the Pyomo Python package.
Pyomo requires an MIP solver, we use the Coin-or Branch and Cut solver (CBC) which can be installed from the command line:

> sudo apt-get install coinor-cbc
