# Dispatch Optimization with Pyomo

We cast the dispatch optimization as a mixed integer linear programming problem and solve it using the Pyomo Python package.
Pyomo requires an MIP solver, we use the Coin-or Branch and Cut solver (CBC) which can be installed from the command line:

> sudo apt-get install coinor-cbc

To test that the installation worked, you can run the `cbc` command from any terminal and see if the solver returns a prompt. 
