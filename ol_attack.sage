from sage.all import *

##############################################################################
# Generate AGCD samples
##############################################################################
def generate_agcd_samples(gamma, eta, rho, n):
    lower_p = 2**(eta - 1)
    upper_p = 2**eta - 1
    p = random_prime(upper_p, lbound=lower_p)
    samples = []
    for _ in range(n):
        r = ZZ.random_element(-2**rho, 2**rho)
        q_bits = gamma - eta
        q = ZZ.random_element(2**(q_bits - 1), 2**q_bits)
        x = p * q + r
        samples.append(x)
    return p, samples

##############################################################################
# Orthogonal lattice attack
#
# Following the standard AGCD orthogonal lattice construction.
#
#          [ x1     x2 ... xt ]
#          [ 2^rho  0  ...  0 ]
#    B  =  [ 0   2^rho ...  0 ]
#          [ 0 ...   0  2^rho ]
#
# The lattice L is spanned by the *columns* of B.
# We reduce B.T (whose rows span L) with BKZ to find t short vectors vᵢ
# orthogonal to w = (1, r1/2^rho, ..., rt/2^rho).
# Collecting those into V and solving V·w = 0 with w[0]=1 recovers the rᵢ's.
##############################################################################
def orthogonal_lattice_attack(samples, rho, block_size=10):
    t = len(samples)

    # Build the (t+1) x t basis matrix B
    B = Matrix(ZZ, t + 1, t)
    for i in range(t):
        B[0, i]     = samples[i]   # first row: the AGCD samples
        B[i + 1, i] = 2**rho       # diagonal block: 2^rho  (Bug 1 fixed: rho not ρ)

    print("[+] Running BKZ on B^T ...")
    # B.T has shape t x (t+1); its rows span L
    V = B.T.BKZ(block_size=block_size)
    print("[+] Reduced basis obtained.")

    # Bug 2 fixed: recover w from the right kernel of V
    # V has shape t x (t+1); w lives in its right kernel
    K = V.right_kernel().basis()
    if len(K) == 0:
        raise ValueError("Right kernel is empty — BKZ did not find enough short vectors.")

    # Normalise so that the first coordinate is 1
    w = K[0]
    if w[0] == 0:
        raise ValueError("w[0] = 0 — unexpected kernel vector shape.")
    w = w / w[0]

    # Recover noise terms rᵢ = round(wᵢ₊₁ · 2^rho)
    noise = [ZZ(round(w[i + 1] * 2**rho)) for i in range(t)]

    # Recover candidate p
    denoised = [samples[i] - noise[i] for i in range(t)]
    p_candidate = gcd(denoised)

    print("[+] Noise terms recovered:", noise[:4], "...")
    print("[+] Candidate p (gcd of denoised samples):", p_candidate)

    # Bug 3 fixed: return only quantities that are actually computed
    return {
        "basis":        B,
        "reduced_basis": V,
        "w":            w,
        "noise":        noise,
        "p_candidate":  p_candidate,
    }

##############################################################################
# Recover p
##############################################################################
def recover_p(samples, rho, block_size=10):
    res = orthogonal_lattice_attack(samples, rho, block_size=block_size)
    g = res["p_candidate"]
    if g == 0:
        print("[-] Recovery failed: gcd is 0.")
        return None
    return g

##############################################################################
# Test
##############################################################################
def test(gamma=157, eta=128, rho=121, n=25, block_size=25):
    print("[+] Generating AGCD samples...")
    p, samples = generate_agcd_samples(gamma, eta, rho, n)
    print("\nSecret p =", p)

    print("\n[+] Launching orthogonal lattice attack...\n")
    res = orthogonal_lattice_attack(samples, rho, block_size=block_size)

    print("\n[+] Attempting recovery of p...\n")
    g = recover_p(samples, rho, block_size=block_size)
    print("\nRecovered value =", g)

    if g is not None and p % g == 0:
        print("\n[SUCCESS] Recovered a divisor of p (likely p itself).")
    else:
        print("\n[FAILURE]")

test()