from math import floor, lgamma, log, log2, pi

from estimator.agcd_parameters import AGCDParameters as Parameters
from estimator import reduction
from estimator.reduction import RC


# Conditions on log2(delta_0) under which an orthogonal lattice (OL) attack with n
# samples works (see the thesis, subsection "Condicoes concretas de sucesso"):
#
#   "typical":   2nd OL attack (lattice L(alpha)), typical sizes:
#                log delta_0 < (eta - rho - (gamma - rho)/n) / n.
#   "pereira":   1st OL attack (lattice orthogonal to x), typical sizes, as in
#                Pereira (PKC 2021): log delta_0 < (eta - rho - gamma/n) / n.
#   "xu_full":   Xu, Sarkar and Hu (2018), condition (7): worst-case bounds under
#                the GSA, all n-1 relations from a single reduction.
#   "xu_single": Xu et al., condition (9) (Remark 2): one relation per reduction,
#                which is repeated 2n-1 times.
#
# The typical-size conditions approximate the number of samples from which the
# attack starts to work; they give lambda. Xu's conditions are sufficient ones:
# they guarantee success, so their cost is only reported as an upper reference,
# and a non-positive Xu bound does not mean that the attack fails.
OL_CONDITIONS = ("typical", "pereira", "xu_full", "xu_single")
OL_FOR_LAMBDA = ("typical", "pereira")
OL_NAMES = {
    "typical": "OL, 2nd attack (typical sizes)",
    "pereira": "OL, 1st attack (Pereira)",
    "xu_full": "OL, Xu et al. (7)",
    "xu_single": "OL, Xu et al. (9)",
}

# Up to this dimension the attack may use HKZ (BKZ with block size n), whose
# root-Hermite factor is taken from the Gaussian heuristic; reduction.beta does
# not return block sizes below 40.
HKZ_MAX_DIM = 40


def _ol_log_delta0(gap, signal, gamma, n, condition):
    """
    Upper bound on log2(delta_0) under which the OL attack works with n samples.

    :param gap: eta - rho.
    :param signal: gamma - rho.
    :param gamma: gamma.
    :param n: number of AGCD samples.
    :param condition: one of ``OL_CONDITIONS``.
    """
    if condition == "typical":
        return (gap - signal / n) / n
    if condition == "pereira":
        return (gap - gamma / n) / n
    if condition == "xu_full":
        return (gap - signal / n - 0.5 * log2(n * (n + 2))) / n
    if condition == "xu_single":
        return (gap - signal / n - 0.5 * log2(4 * n)) / n
    raise ValueError(f"unknown OL condition: {condition}")


def _log2_delta_hkz(n):
    """log2 of the root-Hermite factor of HKZ in dimension n (Gaussian heuristic)."""
    return ((lgamma(n / 2 + 1) / n) / log(2) - 0.5 * log2(pi)) / n


def _entry_bits(gamma, rho, condition):
    """
    Bit size B of the largest basis entries, passed to the cost model (which then
    adds the cost of LLL on entries of that size).

    The 2nd OL attack can use the rounding variant of Xu et al., whose first column
    is floor(x_i / alpha), with alpha ~ 2^rho, so its entries have about gamma - rho
    bits. The basis of the lattice orthogonal to x (1st attack, "pereira") is
    computed from the x_i themselves, which have gamma bits.
    """
    return gamma if condition == "pereira" else gamma - rho


def _log2_cost(red_cost_model, beta, n, condition, B):
    """
    log2 of the cost of the attack: one reduction, or 2n-1 of them for xu_single.

    :param B: bit size of the basis entries (see ``_entry_bits``).
    """
    try:
        log2_T = float(log2(red_cost_model(beta, n, B=B)))
    except (OverflowError, ValueError, ZeroDivisionError):
        return float("inf")
    if condition == "xu_single":
        log2_T += log2(2 * n - 1)
    return log2_T


def _ol_best(gap, signal, gamma, red_cost_model, condition, beta_slack=20):
    """
    Cheapest OL attack under the given condition, minimising over the number of
    samples n and the reduction (HKZ in small dimension, BKZ-beta otherwise).

    The bound on log2(delta_0) rises and then falls with n. The smallest block size
    is reached at the peak, but a slightly larger beta may allow a much smaller n, so
    for each beta in [beta_min, beta_min + beta_slack] we take the smallest n (with
    n >= beta) that reaches it, and keep the cheapest combination.

    :returns: a dict with the attack parameters, or ``None`` if no n works.
    """
    bound = lambda n: _ol_log_delta0(gap, signal, gamma, n, condition)  # noqa: E731
    n_min = max(2, floor(signal / gap) + 1)

    # Walk up the rising side of the bound until it peaks. The bound is N(n)/n,
    # where N(n) rises and then falls (or keeps rising, for the typical-size
    # conditions). If N peaks without becoming positive, the bound is non-positive
    # for every n.
    n_peak, ld_peak = n_min, bound(n_min)
    while True:
        ld = bound(n_peak + 1)
        if ld_peak > 0 and ld < ld_peak:
            break
        if ld_peak <= 0 and ld * (n_peak + 1) < ld_peak * n_peak:
            return None
        n_peak, ld_peak = n_peak + 1, ld

    def smallest_n(target):
        """Smallest n on the rising side with bound(n) >= target (binary search)."""
        lo, hi = n_min, n_peak
        while lo < hi:
            mid = (lo + hi) // 2
            if bound(mid) >= target:
                hi = mid
            else:
                lo = mid + 1
        return lo

    candidates = []  # (n, beta, algorithm)

    # HKZ in small dimension. Every valid n is costed, also past the peak of the
    # bound: the cost models are not monotonic in n at these block sizes.
    for n in range(n_min, HKZ_MAX_DIM + 1):
        if bound(n) >= _log2_delta_hkz(n):
            candidates.append((n, n, "HKZ"))

    # BKZ-beta, which needs a lattice of dimension n >= beta.
    beta_min = reduction.beta(2**ld_peak)
    if beta_min is not None:
        for beta in range(beta_min, beta_min + beta_slack + 1):
            target = log2(float(reduction.delta(beta)))
            if target > ld_peak:
                continue
            # n >= beta may lie past the peak, where the bound decreases again.
            n = max(smallest_n(target), beta)
            if bound(n) < target:
                continue
            candidates.append((n, beta, "BKZ"))

    B = _entry_bits(gamma, gamma - signal, condition)  # rho = gamma - signal
    best = None
    for n, beta, algorithm in candidates:
        log2_T = _log2_cost(red_cost_model, beta, n, condition, B)
        if log2_T == float("inf"):
            continue
        if best is None or log2_T < best["lambda"]:
            best = {
                "condition": condition,
                "n": n,
                "log2_delta0_bound": bound(n),
                "algorithm": algorithm,
                "beta": beta,
                "B": B,
                "T_lattice": 2**log2_T,
                "lambda": log2_T,
            }
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
        Estimate the cost of the known attacks on an AGCD instance.

        lambda is the cheapest of: the OL attacks under the typical-size conditions
        ("typical" and "pereira"), the GCD attack and the trivial attack. The OL
        costs under Xu et al.'s conditions are reported as an upper reference.

        :param params: AGCD parameters.
        :param red_cost_model: How to cost lattice reduction.
        :param jobs: Unused; kept for compatibility with the other estimators.
        :param catch_exceptions: Unused; kept for compatibility.
        :param verbose: Print the intermediate results.
        """
        params = params.normalize()
        gamma, eta, rho = params.gamma, params.eta, params.rho

        # ------------------------------------------------------------------
        # Step 1 - parameter check
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

        if verbose:
            print("\n=== Step 1: Parameters ===")
            print(f"  gamma - rho = {signal:.2f}")
            print(f"  eta   - rho = {gap:.2f}")
            print(f"  Minimum n for the OL attacks: n > {signal / gap:.4f}")

        # ------------------------------------------------------------------
        # Step 2 - OL attacks under each condition, n optimised
        # ------------------------------------------------------------------
        ol = {c: _ol_best(gap, signal, gamma, red_cost_model, c) for c in OL_CONDITIONS}

        if verbose:
            print("\n=== Step 2: OL attacks (n optimised) ===")
            for c, r in ol.items():
                role = "" if c in OL_FOR_LAMBDA else "  [upper reference]"
                if r is None:
                    print(f"  {OL_NAMES[c] + ':':32} no n satisfies the condition{role}")
                else:
                    print(f"  {OL_NAMES[c] + ':':32} n = {r['n']}, {r['algorithm']}"
                          f"-{r['beta']}, log2(T) = {r['lambda']:.2f}{role}")

        # ------------------------------------------------------------------
        # Step 3 - GCD and trivial attacks
        # ------------------------------------------------------------------
        # GCD attack, as per section 4 of "Efficient AGCD-based homomorphic
        # encryption for matrix and vector arithmetic". Vanilla AGCD takes n = 1.
        T_gcd = rho**2 * 2**(rho + rho/2) * gamma * log2(gamma)

        # ------------------------------------------------------------------
        # Step 4 - lambda: the cheapest attack
        # ------------------------------------------------------------------
        inf = float("inf")
        attack_costs = {OL_NAMES[c]: (ol[c]["lambda"] if ol[c] else inf) for c in OL_FOR_LAMBDA}
        best_lattice = min(attack_costs.values())
        attack_costs["GCD"] = log2(T_gcd)
        attack_costs["trivial"] = eta
        best_attack = min(attack_costs, key=attack_costs.get)
        log2_T = attack_costs[best_attack]
        T = 2**log2_T
        lam = log2_T

        feasible = [ol[c] for c in OL_FOR_LAMBDA if ol[c] is not None]
        ol_best = min(feasible, key=lambda r: r["lambda"]) if feasible else None

        if verbose:
            print("\n=== Step 4: Attack costs ===")
            for name, cost in attack_costs.items():
                shown = "n/a" if cost == inf else f"{cost:.2f}"
                print(f"  {name + ':':32} log2(T) = {shown}")
            print(f"  Bottleneck: {best_attack}, log2(T) = {log2_T:.2f} bits")

            print("\n=== Result ===")
            print(f"  T      ~ 2^{log2_T:.2f} clock cycles")
            print(f"  lambda ~ {lam:.2f} bits (set by {best_attack})")
            if eta < best_lattice:
                print("  NOTE: the trivial attack is cheaper than every lattice attack.")
            if log2(T_gcd) < best_lattice:
                print("  NOTE: the GCD attack is cheaper than every lattice attack.")

        return {
            "trivial": {
                "T_trivial (2^eta)": 2**eta,
                "lambda (eta)": eta,
            },
            "gcd": {
                "T_gcd": T_gcd,
                "lambda": log2(T_gcd),
            },
            "ol": ol,
            "ol_bkz": ol_best,
            "final": {
                "lambda": log2_T,
                "T": T,
                "attack": best_attack,
            },
        }


estimate = Estimate()
