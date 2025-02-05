# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sage.all import oo, binomial, log, sqrt, ceil

from .nd import NoiseDistribution, DiscreteGaussian
from .errors import InsufficientSamplesError


@dataclass
class AGCDParameters:
    """The parameters for a Approximated GCD problem instance, ai := p*qi + ri"""
    
    gamma: int  #: the size in bits of ai
    eta: int    #: the size in bits of p
    rho: int    #: the size in bits of ri

    lamda: int  #: expected the security level
    
    Xr: NoiseDistribution = DiscreteGaussian(3.0)  #: the distribution from which the error term ri is drawn

    #: the number of samples allowed to an attacker,
    #: optionally `sage.all.oo` for allowing infinitely many samples.
    m: int = oo

    tag: str = None  #: a name for the parameter set

    def normalize(self):
        # nothing to do
        return self

    def updated(self, **kwds):
        """
        Return a new set of parameters updated according to ``kwds``.
        """
        d = dict(self.__dict__)
        d.update(kwds)
        return AGCDParameters(**d)

 
    def __hash__(self):
        return hash((self.n, self.q, self.Xs, self.Xe, self.m, self.tag))
