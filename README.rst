Security Estimates for AGCD
=======================================

This is the AGCD branch for the celebrated Lattice Estimator.
For general information about the Lattice Estimator, see the `README
for the main branch
<https://github.com/malb/lattice-estimator/blob/main/README.rst>`__.


The main purpose of this estimator is to give designers an easy way to choose parameters resisting known attacks and to enable cryptanalysts to compare their results and ideas with other techniques known in the literature.


Usage examples:
---------------

For using sage notebook (such as `agcd_test_notebook`), run `sage -n
jupyterlab` at the root of the project.

Following is an example for running the estimator for a specified set
of parameters.

  .. code-block:: python
		  
    >>> from estimator import AGCD
    >>> agcd_params = AGCD.Parameters(gamma = 20, eta = 15, rho = 10, lamda = 20)
    >>> AGCD.estimate(agcd_params)

=== Step 1: Feasibility ===
  gamma - rho = 10.00
  eta   - rho = 5.00
  Minimum n for attack: n > 2.0000

=== Step 2: Optimal n and delta_0 ===
  Optimal n (continuous) = 4.0000
  Rounded n              = 4
  log2(delta_0) bound    = 0.051880
  delta_0 bound          = 1.03661465

=== Step 3: BKZ Blocksize ===
  Minimum beta = 40
  delta_0 achieved at beta: 1.01295000

=== Step 4: Attack Cost (Section 3.3, b_max = gamma - rho = 10.00) ===
  Sieving cost:  log2(T_sieve) = 33.08 bits
  LLL cost:      log2(T_lll)   = 2.82 bits
  Bottleneck:    log2(T)        = 33.08 bits

=== Step 5: Asymptotic Cross-check ===
  (gamma-rho)/(eta-rho)^2 * log(...) = -0.5288
  (Should be proportional to lambda = 33.08)

=== Result ===
  T      ≈ 2^33.08 clock cycles
  lambda ≈ 33.08 bits
Time cost for achieving 33.08 bits of security, with beta = 40: 9079715830.98625

{'rho_eff': 10,
 'n': 4,
 'delta0': 1.0366146496280775,
 'beta': 40,
 'T_sieve': 9079715830.98625,
 'T_lll': 7.052685103294118,
 'T': 9079715830.98625,
 'lambda': 33.08,
 'broken': False}

  
         
Evolution
---------

This code is evolving, new results are added and bugs are fixed. Hence, estimations from earlier
versions might not match current estimations. This is annoying but unavoidable. We recommend to also
state the commit that was used when referencing this project.

.. warning :: We give no API/interface stability guarantees. We try to be mindful but we may reorganize the code without advance warning.
