import scipy.optimize as sp
import numpy as np

from math import floor, log2
from estimator.agcd_parameters import AGCDParameters as Parameters
from estimator import reduction
from estimator.reduction import RC


# Variants of the orthogonal lattice (OL) attack of Xu, Sarkar and Hu (2018):
#   "full":   condition (7), the n-1 relations come from a single reduction (i = n-1);
#   "single": condition (9), one relation per reduction (i = 1), repeated over
#             N = 2n samples until N-1 relations are found (Remark 2), so the
#             reduction cost is paid 2n-1 times.
OL_VARIANTS = ("full", "single")


def _ol_log_delta0(gap, signal, n, variant):
    """
    Upper bound on log2(delta_0) under which the OL attack succeeds with n samples.

    :param gap: eta - rho.
    :param signal: gamma - rho.
    :param n: number of AGCD samples (lattice dimension).
    :param variant: "full" (condition (7)) or "single" (condition (9)).
    """
    if variant == "full":
        log_norm = 0.5 * log2(n * (n + 2))
    else:
        log_norm = 0.5 * log2(4 * n)
    return (gap - signal / n - log_norm) / n


def _ol_best(gap, signal, red_cost_model, variant, beta_slack=20):
    """
    Cheapest OL attack of the given variant, minimising over the number of samples n.

    The bound on log2(delta_0) rises and then falls with n. The smallest block size
    beta is reached at the peak, but a slightly larger beta may allow a much smaller
    n, so for each beta in [beta_min, beta_min + beta_slack] we take the smallest n
    that reaches it and keep the cheapest combination.

    :returns: a dict with the attack parameters, or ``None`` if the bound is
              non-positive for every n (the attack does not apply).
    """
    n_min = max(2, floor(signal / gap) + 1)

    # Walk up the rising side of the bound until it peaks. The bound is N(n)/n,
    # where N(n) = gap - signal/n - log_norm(n) rises and then falls. If N peaks
    # without becoming positive, the bound is non-positive for every n (it only
    # creeps towards 0 from below, so waiting for it to fall would never end).
    n_peak = n_min
    ld_peak = _ol_log_delta0(gap, signal, n_peak, variant)
    while True:
        ld = _ol_log_delta0(gap, signal, n_peak + 1, variant)
        if ld_peak > 0 and ld < ld_peak:
            break
        if ld_peak <= 0 and ld * (n_peak + 1) < ld_peak * n_peak:
            return None
        n_peak, ld_peak = n_peak + 1, ld

    beta_min = reduction.beta(2**ld_peak)
    if beta_min is None:
        # delta_0 so close to 1 that no block size is found: treat as infeasible.
        return None
    best = None
    for beta in range(beta_min, beta_min + beta_slack + 1):
        target = log2(float(reduction.delta(beta)))
        if target > ld_peak:
            continue
        # Smallest n on the rising side with bound >= target (binary search).
        lo, hi = n_min, n_peak
        while lo < hi:
            mid = (lo + hi) // 2
            if _ol_log_delta0(gap, signal, mid, variant) >= target:
                hi = mid
            else:
                lo = mid + 1
        n = lo
        try:
            log2_T = float(log2(red_cost_model(beta, n)))
        except OverflowError:  # e.g. ELHL26 computes 2**x with Python floats
            continue
        if log2_T == float("inf"):  # cost beyond float range: block size far too large
            continue
        if variant == "single":
            log2_T += log2(2 * n - 1)
        if best is None or log2_T < best["lambda"]:
            best = {
                "variant": variant,
                "n": n,
                "log2_delta0_bound": _ol_log_delta0(gap, signal, n, variant),
                "beta": beta,
                "T_lattice": 2**log2_T,
                "lambda": log2_T,
            }
    # None also when every cost overflowed (block size far beyond any model).
    return best


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

        small_lattice_dim = 40

        trivial = False
        GCD = False
        broken = False
        hkz = False

        log2_T_lattice = float('inf')

        lambdas = []

        # ------------------------------------------------------------------
        # Step 1 - feasibility check
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
        # Step 2 - optimal n and delta_0 bound (Theorem 2 / eq. 5)
        # ------------------------------------------------------------------
        n_opt = 2 * signal / gap     # optimal n = 2*(gamma-rho)/(eta-rho)
        n = max(2, round(n_opt))
        if n < 40:
            if verbose:
                print(f"  n = {n} is too small, so we may use HKZ instead of BKZ")
            hkz = True
        try:
            log2_T_hkz = 0.292*n + 16.4
            T_hkz = 2**(0.292*n + 16.4)
        except OverflowError:
            log2_T_hkz = float("inf")
            T_hkz = float("inf")

        lambdas.append(log2_T_hkz)
        # Theorem 2 of Xu working condition (eq. 5), solved for log delta_0,
        # using i = n-1 (need n-1 independent orthogonal vectors):
        #
        #   log delta_0 < (1/n) * [ (eta-rho) - (gamma-rho)/n - log sqrt(n*(n-1+3)) ]
        #               = (1/n) * [ gap - signal/n - 0.5*log(n*(n+2)) ]
        #
        #  condition uses i+3 with i = n-1, so i+3 = n+2.
        log_delta0_bound = (1 / n) * (gap - signal/n - 0.5 * log2(n * (n + 2)))
        delta0 = 2 ** log_delta0_bound


        if log_delta0_bound <= 0:
            # Xu's attack cannot achieve a useful delta_0 with this n
            # (scheme is secure against this attack at these parameters)
            if verbose:
                print(f"\n  WARNING: delta_0 bound is non-positive ({log_delta0_bound:.6f}). Falling back to Hilder's estimate?")
                print("  The scheme appears secure against this OL attack.")


        hilder_n = 2*gamma/gap # estimated from Hilder's thesis
        hilder_delta0_bound = 2**(gap/hilder_n - gamma/hilder_n**2) # (gap - gamma/n)/n
        hilder_beta = None
        hilder_log2_T_lattice = float("inf")
        if hilder_delta0_bound > 1:
            hilder_beta = reduction.beta(hilder_delta0_bound)
            hilder_log2_T_lattice = log2(red_cost_model(hilder_beta, hilder_n))
            lambdas.append(hilder_log2_T_lattice)

        if verbose:
            print("######## Hilder's estimation ########")
            print(f"  hilder_n               = {hilder_n}")
            print(f"  hilder_delta0_bound    = {hilder_delta0_bound:.8f}")
            if hilder_delta0_bound > 1:
                print(f"  hilder beta            = {hilder_beta:.8f}")
                print(f"  hilder lambda          = {hilder_log2_T_lattice:.8f}")
            print("######## ######## ######## ########")

        # ------------------------------------------------------------------
        # Step 3 - OL attack: both variants of Xu et al., optimising n
        # ------------------------------------------------------------------
        ol = {variant: _ol_best(gap, signal, red_cost_model, variant) for variant in OL_VARIANTS}
        feasible = [r for r in ol.values() if r is not None]
        if feasible:
            ol_best = min(feasible, key=lambda r: r["lambda"])
            log2_T_lattice = ol_best["lambda"]
            lambdas.append(log2_T_lattice)
        else:
            ol_best = None

        if verbose:
            print(f"\n=== Step 3: OL attack (Xu et al.), n optimised ===")
            print(f"  Closed-form n (eq. 8)  = {float(n_opt):.4f}")
            for variant, r in ol.items():
                if r is None:
                    print(f"  {variant:>6}: bound non-positive for every n, attack does not apply")
                else:
                    print(f"  {variant:>6}: n = {r['n']}, beta = {r['beta']}, "
                          f"log2(T) = {r['lambda']:.2f}")

        # Estimates time taken for GCD attack (as per section 4 of "Efficient AGCD-based homomorphic encryption for matrix and vector arithmetic)
        # Vanilla AGCD takes n = 1
        T_gcd = rho**2 * 2**(rho + rho/2) * gamma * log2(gamma)
        lambdas.append(log2(T_gcd))
        lambdas.append(eta)

        # Compare the GCD and trivial attacks against the cheapest lattice attack
        # (OL of Xu et al., Hilder's estimate or HKZ), not only against Xu's OL.
        lattice_costs = {
            "OL (Xu et al.)": log2_T_lattice,
            "OL (Hilder)": hilder_log2_T_lattice,
            "HKZ": log2_T_hkz,
        }
        best_lattice = min(lattice_costs.values())
        GCD = log2(T_gcd) < best_lattice
        trivial = eta < best_lattice
        attack_costs = dict(lattice_costs, **{"GCD": log2(T_gcd), "trivial": eta})
        best_attack = min(attack_costs, key=attack_costs.get)

        # Actual running time of attack is the minimum of the attacks taken into consideration
        #log2_T = min(log2_T_lattice, eta, log2(T_gcd)) if delta0 > 1 else  min(eta, log2(T_gcd))
        log2_T = min(lambdas)
        T = 2**log2_T

        # lambda in bits
        lam = log2_T

        if verbose:
            print(f"\n=== Step 4: Attack Cost  ===")
            for name, cost in attack_costs.items():
                shown = "n/a" if cost == float("inf") else f"{cost:.2f}"
                print(f"  {name + ':':16} log2(T) = {shown}")
            print(f"  Bottleneck:      {best_attack}, log2(T) = {log2_T:.2f} bits")

            print(f"\n=== Result ===")
            print(f"  T      ~ 2^{log2_T:.2f} clock cycles")
            print(f"  lambda ~ {lam:.2f} bits (set by {best_attack})")
            if trivial:
                print(f"  NOTE: the trivial attack is cheaper than every lattice attack.")
            if GCD:
                print(f"  NOTE: the GCD attack is cheaper than every lattice attack.")

        return {"trivial": {
                            "T_trivial (2^eta)": 2**eta,
                            "lambda (eta)": eta
                           },
                "gcd": {
                        "T_gcd": T_gcd,
                        "lambda": log2(T_gcd)
                       },
                "ol_bkz": ol_best,
                "ol_full": ol["full"],
                "ol_single": ol["single"],
                "ol_hilder_fallback": {
                    "hilder_n": hilder_n,
                    "hilder_delta0_bound": hilder_delta0_bound,
                    "hilder_beta": hilder_beta,
                    "lambda": hilder_log2_T_lattice,
                },
                "ol_hkz": {"T": T_hkz, "lambda": log2_T_hkz},
                "final": {
                            "lambda": log2_T,
                            "T": T
                         }
               }

estimate = Estimate()
