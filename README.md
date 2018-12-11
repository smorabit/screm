## Stochastic bioChemical REaction network simulations with continuous Markov processes

screm.py is a Python module containing functions suitable for modeling biochemical reaction networks using stochastic simulations. This was created as a final project for Math 227A (Mathematical Computational Biology) at UC Irvine. This framework has been applied to the glycolysis reaction cascade, which can be seen in the `glycolysis.ipynb` notebook in this repo.

### Installation

1. Intall anaconda for Python3. It should contain all required dependencies.
2. Clone this repo. Import screm.py using the following syntax:
```python
from screm import Screm
scr = Screm()
```

### Functions

* **delta_t(l)**
    * Computes transition time using the equation $ \Delta t = -\frac{\ln\sigma}{\lambda} $ where $ \sigma $ is a uniform random variable, and $ \lambda $ is the reaction coefficient, which is computed as $ \lambda = k_i \Pi_{\rm j=1}^{\rm n} [N_j] $ for reaction $i$ where $k_i$ is the rate of that reaction and $[N_j]$ is the number of molecules of reactant j for all $(j=1\rightarrow n)$ in reaction $i$.
    
    
   
* **simulate_reactions(v, reaction_net, iterations, verbose=False)**
    * Simulates stochastic chemical reactions within a defined reaction network, given initial molecule counts vector $v$ and a defined number of iterations using a continuous markov process. 
    * Output:
        * T: array of each transition time
        * V: A matrix containing the number of molecules for each reactant for each transition in the markov process.
        * I: final number of iterations
    
    
* **plot_simulation(T,V,I, reaction_net, metabolites="all", figsize=(7,3), dpi=170)**
    * Plots the number of molecules versus in silico time for all specified metabolites in a defined reaction network.
   
   
* **plot_simulations_violin(model_results, init, ax, datatype="iter", poly=1)**
    * For a set of repeated simulations with the same conditions, plot violin plots of time to simulation convergence (or number of iterations to convergence) as well as a best fit line.
    
    
* **plot_residuals(model_results, init, ax, datatype="iter", poly=1)**
    * Plots the difference between the best fit polynomail line and the actual data points.

### Glycolysis example

We aim to investigate the glycolysis reaction cascade using a stochastic in silico framework. For all subsequent models, we begin with constructing a reaction network using the following steps:

1. Create a list of reactants to be modeled in the reaction network.
2. Construct a stoichiometry matrix for the reaction, in which each row is a reactant, and each column is a reaction.
3. Create a dictionary where the key is the reaction number, and the value is a list of reactants that that reaction depends on. For example, the first reaction in glycolysis depends on presence of Glucose, so ```reaction_relation[0] = ["Glucose"]```
4. Perform a literature search to find reaction rates for all relecant reactions.
5. Compile these features into a reaction_net dictionary object.

Please look at the `glycolysis.ipynb` for the rest of the analysis.
