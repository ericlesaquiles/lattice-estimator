import scipy.optimize as sp
import numpy as np

from math import log2
from estimator.agcd_parameters import AGCDParameters as Parameters
from estimator import reduction
from estimator.reduction import RC


class Estimate:

    def __call__(
        self,
        params,
        red_cost_model=RC.MATZOV,
        jobs=1,
        catch_exceptions=True,
        verbose=True,
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
        gamma, eta, rho = params.gamma, params.eta, params.rho

        trivial = False
        GCD = False
        broken = False

        log2_T_bkz = float('inf')

        # ------------------------------------------------------------------
        # Step 1 — feasibility check
        # ------------------------------------------------------------------
        if eta <= rho:
            raise ValueError(
                f"Invalid parameters: eta ({eta}) must be greater than rho ({rho}). "
                "The attack cannot work."
            )
        if gamma <= eta:
            raise ValueError(
                f"Invalid parameters: gamma ({gamma}) must be greater than eta ({eta})."
            )

        gap = eta - rho          # eta - rho
        signal = gamma - rho     # gamma - rho

        min_n = signal / gap         # minimum n for attack to be possible

        if verbose:
            print("\n=== Step 1: Feasibility ===")
            print(f"  gamma - rho = {signal:.2f}")
            print(f"  eta   - rho = {gap:.2f}")
            print(f"  Minimum n for attack: n > {float(min_n):.4f}")

        # ------------------------------------------------------------------
        # Step 2 — optimal n and delta_0 bound (Theorem 2 / eq. 5)
        # ------------------------------------------------------------------
        n_opt = 2 * signal / gap     # optimal n = 2*(gamma-rho)/(eta-rho)
        n = max(2, round(n_opt))

        # Theorem 2 of Xu working condition (eq. 5), solved for log delta_0,
        # using i = n-1 (need n-1 independent orthogonal vectors):
        #
        #   log delta_0 < (1/n) * [ (eta-rho) - (gamma-rho)/n - log sqrt(n*(n-1+3)) ]
        #               = (1/n) * [ gap - signal/n - 0.5*log(n*(n+2)) ]
        #
        #  condition uses i+3 with i = n-1, so i+3 = n+2.
        log_delta0_bound = (1 / n) * (gap - signal/n - 0.5 * log2(n * (n + 2)))

        if log_delta0_bound <= 0:
            # Xu's attack cannot achieve a useful delta_0 with this n
            # (scheme is secure against this attack at these parameters)
            if verbose:
                print(f"\n  WARNING: delta_0 bound is non-positive ({log_delta0_bound:.6f}). Falling back to Hilder's estimate?")
                print("  The scheme appears secure against this OL attack.")



        hilder_n = 2*gamma/gap # estimated from Hilder's thesis
        hilder_delta0_bound = 2**(gap/n - gamma/n**2) # (gap - gamma/n)/n
        if verbose:
            print("######## Hilder's estimation ########")
            print(f"  hilder_n               = {hilder_n}")
            print(f"  hilder_delta0_bound    = {hilder_delta0_bound:.8f}")
            if hilder_delta0_bound > 1:
                beta = reduction.beta(hilder_delta0_bound)
                hilder_log2_T_bkz = log2(red_cost_model(beta, hilder_n))
                print(f"  hilder beta            = {beta:.8f}")
                print(f"  hilder lambda          = {hilder_log2_T_bkz:.8f}")
            print("######## ######## ######## ########")




        delta0 = 2 ** log_delta0_bound

        if delta0 > 1:
            if verbose:
                print(f"\n=== Step 2: Optimal n and delta_0 ===")
                print(f"  Optimal n (continuous) = {float(n_opt):.4f}")
                print(f"  Rounded n              = {n}")
                print(f"  log2(delta_0) bound    = {log_delta0_bound:.6f}")
                print(f"  delta_0 bound          = {delta0:.8f}")
                print(f"  hilder_n               = {hilder_n}")
                print(f"  hilder_delta0_bound    = {hilder_delta0_bound:.8f}")

            # ------------------------------------------------------------------
            # Step 3 — minimum BKZ blocksize beta
            # ------------------------------------------------------------------

            beta = reduction.beta(delta0)
            if beta is None:
                raise RuntimeError(
                    "Could not find a valid beta. The scheme may be extremely secure "
                    "or the parameters are unusual."
                )

            if verbose:
                print(f"\n=== Step 3: BKZ Blocksize ===")
                print(f"  Minimum beta = {beta}")
                print(f"  delta_0 achieved at beta: {reduction.delta(beta):.8f}")

            # Sieving-dominated cost: T = 8 * n * 2^(0.292*beta + 16.4)
            # log2_T_sieve = math.log2(8 * n) + 0.292 * beta + 16.4
            # Actually we outsource the cost model
            log2_T_bkz = log2(red_cost_model(beta, n))

            # Actual cost is the bottleneck
            log2_T = log2_T_bkz

        # Estimates time taken for GCD attack (as per section 4 of "Efficient AGCD-based homomorphic encryption for matrix and vector arithmetic)
        # Vanilla AGCD takes n = 1
        T_gcd = rho**2 * 2**(rho + rho/2) * gamma * log2(gamma)

        # Return if GCD attack is better than reduction
        if delta0 < 1 or log2(T_gcd) < log2_T:
            GCD = True
        if delta0 < 1 or eta < log2_T:
            trivial = True

        # Actual running time of attack is the minimum of the attacks taken into consideration
        log2_T = min(log2_T, eta, log2(T_gcd)) if delta0 > 1 else  min(eta, log2(T_gcd))
        T = 2**log2_T

        # lambda in bits
        lam = log2_T

        if verbose:
            print(f"\n=== Step 4: Attack Cost  ===")
            if delta0 > 1:
                print(f"  Reduction cost:  log2(T_bkz) = {log2_T_bkz:.2f} bits")
            print(f"  Trivial cost:  eta           = {eta} bits")
            print(f"  GCD cost:      log2(T_gcd)   = {log2(T_gcd):.2f} bits")

            print(f"  Bottleneck:    log2(T)       = {log2_T:.2f} bits")

            print(f"\n=== Result ===")
            print(f"  T      ≈ 2^{log2_T:.2f} clock cycles")
            print(f"  lambda ≈ {lam:.2f} bits")
            if trivial:
                print(f"  WARNING: Trivial attack is better than lattice reduction.")
            if GCD:
                print(f"  WARNING: GCD attack is better than lattice reduction.")

        return {
            "delta0":   delta0,
            "beta":     beta if delta0 > 1 else None,
            "T_bkz":  2**log2_T_bkz,
            "T":        2**log2_T,
            "lambda":   lam,
        }

estimate = Estimate()
