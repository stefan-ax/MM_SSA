# MM_SSA
Bayesian Inference for Michaelis-Menten Enzyme Kinetics

## How to use?
Go for mm_ssa.py.

```
class MM_SSA():
    '''
    Gillespie SSA implementation for Michaelis Menten enzyme kinetics
    
    E + S -> C (k1)
    C -> E + S (k2)
    C -> E + P (k3)
    
    Input parameters:
        E0
        S0
        C0
        P0
        k1
        k2
        k3
        max_time (default set to 100)
        initial_time (default set to 0)
        
    Default usage:
        model=MM_SSA(1, 10, 1, max_time =20, E0 = 100, S0 = 10, C0 = 0, P0 = 0)
        model.multi_plot(reactant = 'E')
    '''
    
    __author__ = 'Stefan Alexandru Obada'
    __date__ = '09/06/2019'
    
```
