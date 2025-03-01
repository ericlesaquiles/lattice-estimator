import scipy.optimize as sp
import numpy as np

from estimator.agcd_parameters import AGCDParameters as Parameters
from estimator import reduction

from .conf import (
    red_cost_model as red_cost_model_default,
    red_shape_model as red_shape_model_default,
)


class Estimate:

    def _max_delta_estimate(self, eta, gamma, rho):
        er = eta - rho
        gr = gamma - rho
        return 2**( er**2/(4*gr) - er/(2*gr) - er/(4*gr)*np.log((gr/er)**2 + gr/er))

    def _estimate_beta_from_delta(delta):
        beta = sp.newton(lambda n: ((beta/(2*np.pi*np.e)) * (np.pi*np.e)**(1/beta))**(1/(2*(beta - 1))), 1)
        return beta
    
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

        delta0 = self._max_delta_estimate(eta, gamma, rho)
        beta = reduction.beta(delta0)
        n = sp.newton(lambda n: (gamma - rho)/n - (eta - rho) + n*np.log(delta0) + np.log(np.sqrt(n**2 +2*n)), 1) + 1
        
        # TODO move this block to orthogonal_agcd.py
        bkz_cost_constant = 8 # See justification for this number at page 10 of https://eprint.iacr.org/2017/047.pdf

        cost = bkz_cost_constant * n * 2**(0.292*beta + 16.4)

        print(f'Cost for beta = {beta}: {cost}, with delta: {delta0!r}')
            
estimate = Estimate()
