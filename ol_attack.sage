from sage.all import *

##############################################################################
# Generate AGCD samples
##############################################################################

gamma = 157
eta   = 128
rho   = 121

def generate_agcd_samples(gamma, eta, rho, n):

    lower_p = 2^(eta - 1)
    upper_p = 2^eta - 1

    p = random_prime(upper_p, lbound=lower_p)

    samples = []

    for _ in range(n):

        r = ZZ.random_element(-2^rho, 2^rho)

        q_bits = gamma - eta
        q = ZZ.random_element(2^(q_bits - 1), 2^q_bits)

        x = p*q + r

        samples.append(x)

    return p, samples


##############################################################################
# Orthogonal lattice attack
#
# Following the standard AGCD orthogonal lattice construction.
##############################################################################

def orthogonal_lattice_attack(samples, block_size = 10):

    t = len(samples)

    # We use:
    #
    #          [ 1  0  ...  0  x1 ]
    #          [ 0  1  ...  0  x2 ]
    #    B  =  [            ...   ]
    #          [ 0  0  ...  1  xt ]
    #          [ 0  0  ...  0  X  ]
    #
    # where X is a scaling factor.
    #
    # Short vectors correspond to relations:
    #
    #    sum(a_i x_i) ≈ small
    #

    X = 2^(samples[0].nbits())

    B = Matrix(ZZ, t + 1, t + 1)

    for i in range(t):
        B[i, i] = 2**rho
        B[i, t] = samples[i]

    B[t, t] = X

    print("[+] Running BKZ...")

    Bred = B.BKZ(block_size = block_size)

    print("[+] Reduced basis obtained.")

    short = Bred[0]

    coeffs = list(short[:t])

    # Compute the integer combination
    relation = sum(coeffs[i] * samples[i] for i in range(t))

    print("\n[+] Short relation found:")
    print("coefficients =", coeffs)

    print("\n[+] Combination:")
    print("R =", relation)

    return {
        "basis": B,
        "reduced_basis": Bred,
        "coefficients": coeffs,
        "relation": relation
    }


##############################################################################
# Recover p using several short relations
##############################################################################

def recover_p(samples, trials=10):

    relations = []

    for _ in range(trials):

        res = orthogonal_lattice_attack(samples)

        R = abs(res["relation"])

        if R != 0:
            relations.append(R)

    if len(relations) < 2:
        print("[-] Not enough relations.")
        return None

    g = relations[0]

    for r in relations[1:]:
        g = gcd(g, r)

    return g


##############################################################################
# Test
##############################################################################


def test(gamma = 157,
         eta   = 128,
         rho   = 121,
         n     = 9):
    print("[+] Generating AGCD samples...")
    
    p, samples = generate_agcd_samples(gamma, eta, rho, n)
    
    print("\nSecret p =", p)
    
    print("\n[+] Launching orthogonal lattice attack...\n")
    
    res = orthogonal_lattice_attack(samples)
    
    print("\n[+] Attempting recovery of p...\n")
    
    g = recover_p(samples, trials=5)
    
    print("\nRecovered gcd =", g)
    
    if g % p == 0:
        print("\n[SUCCESS] Recovered a multiple of p.")
    else:
        print("\n[FAILURE]")