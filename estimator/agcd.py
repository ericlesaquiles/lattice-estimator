import scipy.optimize as sp
import numpy as np

from estimator.agcd_parameters import AGCDParameters as Parameters
from estimator.reduction import delta

from .conf import (
    red_cost_model as red_cost_model_default,
    red_shape_model as red_shape_model_default,
)


class Estimate:
    
    def __call__(
        self,
        params,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        deny_list=tuple(),
        add_list=tuple(),
        jobs=1,
        catch_exceptions=True,
    ):
        """
        Run all estimates, based on the default cost and shape models for lattice reduction.

        :param params: AGCD parameters.
        :param red_cost_model: How to cost lattice reduction. 
        :param deny_list: skip these algorithms
        :param add_list: add these ``(name, function)`` pairs to the list of algorithms to estimate.a
        :param jobs: Use multiple threads in parallel.
        :param catch_exceptions: When an estimate fails, just print a warning.
        """
        params = params.normalize()
        gamma, eta, rho, lamda = params.gamma, params.eta, params.rho, params.lamda
        
        # TODO move this block to orthogonal_agcd.py
        beta_range = range(20, 200, 10) # TODO check whether this testing range is ok
        bkz_cost_constant = 8 # See justification for this number on the paper
        for beta in beta_range:
            delta0 = delta(beta)
            
            n = sp.newton(lambda n: (gamma - rho)/n - (eta - rho) + n*np.log(delta0) + np.log(np.sqrt(n**2 +2*n)), 1) + 1
            cost = bkz_cost_constant * n * 2**(0.292*beta + 16.4)

            print(f'Cost for beta = {beta}: {cost}, with delta: {delta0!r}')
            
estimate = Estimate()
