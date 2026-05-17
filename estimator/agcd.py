import scipy.optimize as sp
import numpy as np
import math

from estimator.agcd_parameters import AGCDParameters as Parameters
from estimator import reduction
from estimator.reduction import RC


class Estimate:

    def __call__(
        self,
        params,
        red_cost_model=RC.ELHL26,
        deny_list=tuple(),
        add_list=tuple(),
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
        rho_eff = rho

        # ------------------------------------------------------------------
        # Step 1 — feasibility check
        # ------------------------------------------------------------------
        if eta <= rho_eff:
            raise ValueError(
                f"Invalid parameters: eta ({eta}) must be greater than rho_eff ({rho_eff}). "
                "The attack cannot work."
            )
        if gamma <= eta:
            raise ValueError(
                f"Invalid parameters: gamma ({gamma}) must be greater than eta ({eta})."
            )

        gap = eta - rho_eff          # eta - rho
        signal = gamma - rho_eff     # gamma - rho

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
        n = max(2, round(n_opt))     # must be an integer >= 2

        # Theorem 2 working condition (eq. 5), solved for log delta_0,
        # using i = n-1 (need n-1 independent orthogonal vectors):
        #
        #   log delta_0 < (1/n) * [ (eta-rho) - (gamma-rho)/n - log sqrt(n*(n-1+3)) ]
        #               = (1/n) * [ gap - signal/n - 0.5*log(n*(n+2)) ]
        #
        # Note: condition uses i+3 with i = n-1, so i+3 = n+2.

        log_delta0_bound = (1 / n) * (gap - signal / n - 0.5 * math.log2(n * (n + 2)))

        if log_delta0_bound <= 0:
            # Attack cannot achieve a useful delta_0 with this n
            # (scheme is secure against this attack at these parameters)
            if verbose:
                print(f"\n  WARNING: delta_0 bound is non-positive ({log_delta0_bound:.6f}).")
                print("  The scheme appears secure against this OL attack.")

        delta0 = 2 ** log_delta0_bound

        if verbose:
            print(f"\n=== Step 2: Optimal n and delta_0 ===")
            print(f"  Optimal n (continuous) = {float(n_opt):.4f}")
            print(f"  Rounded n              = {n}")
            print(f"  log2(delta_0) bound    = {log_delta0_bound:.6f}")
            print(f"  delta_0 bound          = {delta0:.8f}")

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

        # ------------------------------------------------------------------
        # Step 4 — concrete attack cost
        # Section 3.3 uses b_max = gamma - rho_eff (rounding reduces entries)
        # ------------------------------------------------------------------
        b_max = signal   # gamma - rho_eff

        # Sieving-dominated cost: T = 8 * n * 2^(0.292*beta + 16.4)
        # log2_T_sieve = math.log2(8 * n) + 0.292 * beta + 16.4
        # Actually, we might want to outsource the cost model
        log2_T_bkz = red_cost_model(beta, n)

        # LLL-dominated cost: T = 0.00127 * n^3.18 * b_max^1.83
        T_lll = 0.00127 * (n ** 3.18) * (b_max ** 1.83)
        log2_T_lll = math.log2(T_lll) if T_lll > 0 else 0

        # Actual cost is the bottleneck
        log2_T = max(log2_T_bkz, log2_T_lll)


        # Estimates time taken for GCD attack (as per section 4 of "Efficient AGCD-based homomorphic encryption for matrix and vector arithmetic)
        # Vanilla AGCD takes n = 1
        T_gcd = rho**2 * 2**(rho + rho/2) * gamma * math.log2(gamma)

        # Return if GCD attack is better than reduction
        if math.log2(T_gcd) < log2_T:
            GCD = True
        if eta < log2_T:
            trivial = True

        # Actual running time of attack is the minimum of the attacks taken into consideration
        log2_T = min(log2_T, eta, math.log2(T_gcd))

        lam = log2_T   # lambda in bits

        if verbose:
            print(f"\n=== Step 4: Attack Cost  ===")
            print(f"  Sieving cost:  log2(T_sieve) = {log2_T_bkz:.2f} bits")
            print(f"  LLL cost:      log2(T_lll)   = {log2_T_lll:.2f} bits")
            print(f"  Reduction regime dominated by {"sieving" if log2_T_lll < log2_T_bkz else "LLL"}")
            print(f"  Trivial cost:  eta           = {eta} bits")
            print(f"  GCD cost:      log2(T_gcd)   = {math.log2(T_gcd):.2f} bits")

            print(f"  Bottleneck:    log2(T)       = {log2_T:.2f} bits")

            print(f"\n=== Result ===")
            print(f"  T      ≈ 2^{log2_T:.2f} clock cycles")
            print(f"  lambda ≈ {lam:.2f} bits")
            if trivial:
                print(f"  WARNING: Trivial attack is better than lattice reduction.")
            if GCD:
                print(f"  WARNING: GCD attack is better than lattice reduction.")
            if broken:
                print(f"  WARNING: Scheme is BROKEN by implicit factorization check.")



        print(f'Estimated time cost for achieving {lam:.2f} bits of security, with beta = {beta}: {2**log2_T} clock cycles')

        return {
            "rho_eff":  rho_eff,
            "n":        n,
            "delta0":   delta0,
            "beta":     beta,
            "T_sieve":  2**log2_T_bkz,
            "T_lll":    T_lll,
            "T":        2**log2_T,
            "lambda":   lam,
            "broken":   broken,
        }

estimate = Estimate()
