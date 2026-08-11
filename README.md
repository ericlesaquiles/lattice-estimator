**Security Estimates for AGCD (Approximate Greatest Common Divisor) — a fork of the Lattice Estimator**

This artifact accompanies the paper *"On the concrete hardness of the Approximate GCD problem "*, submitted to SBSeg 2026, main track. The paper studies the security of cryptographic schemes based on the AGCD (*Approximate Greatest Common Divisor*) problem, focusing on the orthogonal lattice attack by Xu, Sarkar and Hu, and proposes a cost estimator for this attack, in the spirit of the well-established [Lattice Estimator](https://github.com/malb/lattice-estimator) for LWE/SIS/NTRU problems.

The artifact consists of:

- a **security estimator** (`estimator/agcd.py`, `estimator/agcd_parameters.py`) that, given a set of AGCD parameters `(gamma, eta, rho)`, computes the estimated cost (in bits, `log2(T)`) of the orthogonal lattice attack, the trivial attack, and the GCD attack, following Chen's BKZ blocksize analysis and the lattice-reduction cost models already implemented in the Lattice Estimator (MATZOV, ELHL26, among others);
- a **practical implementation of the orthogonal lattice attack** in SageMath (`ol_attack.sage`), used to empirically validate the estimator's predictions on toy instances;
- a **Jupyter/SageMath notebook** (`agcd_test_notebook.ipynb`) applying the estimator to parameter sets taken from the literature (DGHV, the Cheon et al. scheme, Benarroch–Lepoint, among others);
- a parameter-sweep utility (`param_sweep.py`) reused from the original Lattice Estimator.

The goal of the artifact is to let reviewers (i) reproduce the security estimates presented in the paper for the AGCD parameters studied, and (ii) experimentally validate that the implemented orthogonal lattice attack recovers the secret `p` on small instances, consistent with the theoretical analysis.

# Structure of the readme.md

This `README.md` follows the minimal template required by the SBSeg 2026 Artifact Evaluation Committee, containing:

1. **Project Title** — identification of the artifact and the associated paper.
2. **Structure of the readme.md** — this index.
3. **Seals Considered** — quality seals the authors are requesting for this artifact.
4. **Basic information** — execution environment, hardware and software requirements.
5. **Dependencies** — libraries, versions and third-party resources needed.
6. **Security concerns** — risks (or lack thereof) in running the artifact.
7. **Installation** — step-by-step instructions to obtain and prepare the artifact for execution.
8. **Minimal test** — a quick check that the installation worked.
9. **Experiments** — step-by-step instructions to reproduce the paper's claims, organized by claim.
10. **LICENSE** — the license under which the artifact is distributed.

The repository's directory structure is as follows:

```
.
├── estimator/                  # Lattice Estimator source code (fork), including the AGCD module
│   ├── agcd.py                 # Core of the AGCD security estimator (orthogonal lattice attack,
│   │                           #   trivial attack and GCD attack)
│   ├── agcd_parameters.py      # AGCD parameter data structure (gamma, eta, rho, ...)
│   ├── reduction.py            # Lattice-reduction cost models (MATZOV, ELHL26, ...)
│   └── ...                     # Remaining modules inherited from the Lattice Estimator (LWE, SIS, NTRU, etc.)
├── ol_attack.sage              # Practical implementation of the orthogonal lattice attack in SageMath
├── agcd_test_notebook.ipynb    # Notebook applying the estimator to parameters from the literature
├── param_sweep.py              # Utility for parameter sweeps and plot generation
├── docker/
│   ├── Dockerfile              # Image based on the upstream malb/lattice-estimator (do NOT use for this artifact)
│   ├── Dockerfile.dev          # Image that uses this fork's local code (recommended)
│   └── README.md               # Container-specific usage instructions
├── docs/                       # Sphinx documentation (inherited from the Lattice Estimator)
├── requirements.txt            # Additional Python dependencies (lint, tests, documentation)
├── pyproject.toml              # black/flake8/pytest configuration (includes doctests)
├── Makefile                    # Target for building the Sphinx documentation
└── README.rst                  # Original technical README of the fork (AGCD API usage example)
```

# Seals Considered

The seals considered are: **Available, Functional and Sustainable**.

- **Available**: the code is publicly accessible on GitHub, on the [`agcd`](https://github.com/ericlesaquiles/lattice-estimator/tree/agcd) branch of the `ericlesaquiles/lattice-estimator` repository.
- **Functional**: the estimator can be run locally or via a Docker container, and reviewers can observe its functionality through the minimal test and the experiments described below.
- **Sustainable**: the code is modularized following the original Lattice Estimator's architecture (one module per problem/technique), documented with docstrings, and the files relevant to each paper claim are explicitly indicated in the Experiments section.

The **Reproducible Experiments** seal is not being requested in this submission because the paper's full numerical results depend on a parameter sweep that may take several hours on commodity hardware.

# Basic information

This section describes the environment needed to run and replicate the experiments.

**Hardware:**

- CPU: any modern x86-64 processor (the code does not use a GPU).
- RAM: at least 2 GB free for the minimal test and the notebook's example parameters; 4–8 GB is recommended for broader parameter sweeps with `param_sweep.py` or for larger BKZ blocksizes in `ol_attack.sage`.
- Disk: approximately 5 GB free (the `sagemath/sagemath` image is large, ~3–4 GB).
- Time: installing the Docker image with SageMath can take 15 to 40 minutes depending on the internet connection; running individual experiments takes from seconds to a few minutes (see per-claim estimates below).

**Software:**

- Operating system: any system with Docker support (Linux, macOS, Windows with WSL2), **or** Linux/macOS with SageMath installed natively.
- [SageMath](https://www.sagemath.org/) ≥ 9.5 (`Dockerfile.dev` uses the `sagemath/sagemath:9.5` image; a more recent native installation also works).
- Python 3.9+ (bundled with the SageMath distribution).
- Docker or Podman, if the recommended installation route (container) is chosen.
- Internet access is only needed during installation, to download the SageMath base image and Python dependencies via `pip`. No internet access is required while running the experiments.

# Dependencies

Additional Python dependencies (beyond what already ships with SageMath) are listed in [`requirements.txt`](../requirements.txt):

```
black
flake8
pyproject-flake8
pytest
pytest-xdist
nbmake
sphinx-book-theme
sphinxcontrib-jupyter
```

These dependencies are mainly used for linting, automated testing (`pytest`, including doctests and notebook execution via `nbmake`), and documentation generation — they are not strictly required just to run the estimator interactively, but they are required to reproduce the full test suite.

The remaining functionality (lattice algebra, BKZ, cost estimation) relies entirely on the libraries already bundled with SageMath (there is no need to separately install `fpylll`, `NTL`, etc. — they already ship with the `sagemath/sagemath` image).

No external benchmarks, third-party datasets, API keys, or credentials are needed to run any part of this artifact: all parameters used in the experiments (DGHV, the Cheon et al. scheme, Benarroch–Lepoint) are embedded in the `agcd_test_notebook.ipynb` notebook and/or described in the paper, and the AGCD samples used in `ol_attack.sage` are generated synthetically at runtime.

**Versions used during development:**

- SageMath 9.5+ (as configured in `docker/Dockerfile.dev`)
- Python 3.9+ (bundled with SageMath 9.5)

# Security concerns

Running this artifact **poses no risk** to reviewers:

- No network access is required beyond the initial package installation (via `pip`/`apt`, during the Docker image build process).
- No personal data is collected, transmitted, or exposed.
- No credentials, API keys, or access to third-party services are needed.
- All processing is local and bounded by CPU/RAM; the only caution is to avoid running parameter sweeps (`param_sweep.py`) with very wide ranges or very high BKZ blocksizes in `ol_attack.sage`, which can consume significant CPU time and memory (mitigated by using the small, purpose-chosen test parameters suggested in the Experiments section).
- We recommend running the artifact inside the provided Docker container, isolating it from the rest of the reviewer's system.

# Installation

The recommended process uses Docker, avoiding the need to manually install SageMath.

1. Clone the repository, on the `agcd` branch:

   ```bash
   git clone --branch agcd --single-branch https://github.com/ericlesaquiles/lattice-estimator.git
   cd lattice-estimator
   ```

2. Build the Docker image from the local code (uses `docker/Dockerfile.dev`, which embeds this fork's code — **do not** use `docker/Dockerfile`, which clones the upstream `malb/lattice-estimator` repository and does not contain the AGCD module):

   ```bash
   docker build -t agcd-estimator -f docker/Dockerfile.dev .
   ```

   This step downloads the `sagemath/sagemath:9.5` base image (~3–4 GB) and installs the dependencies from `requirements.txt`. Expected time: 15–40 minutes, depending on the connection.

3. Start an interactive container with SageMath:

   ```bash
   docker run -it --name agcd-estimator-container agcd-estimator sage
   ```

   At the end of this step, you should be in an interactive SageMath shell (`sage:`) with the `/lattice-estimator` directory available and the `estimator` package importable.

**Alternative (native installation, without Docker):** install [SageMath ≥ 9.5](https://doc.sagemath.org/html/en/installation/index.html) for your operating system, then, from the root of the cloned repository, run `sage -pip install -r requirements.txt`. All commands below assume you have a `sage` shell available (native or inside the container).

# Minimal test

This minimal test confirms that the installation worked and lets reviewers observe the artifact's core functionality: estimating the security of an AGCD instance.

1. Inside the `sage` shell (native or in the container, with the working directory at the repository root), run:

   ```python
   from estimator import AGCD
   agcd_params = AGCD.Parameters(gamma=20, eta=15, rho=10)
   AGCD.estimate(agcd_params)
   ```

2. **Expected result:** the command prints a step-by-step report (attack feasibility, minimum BKZ blocksize, estimated sieving/LLL cost) and returns a Python dictionary containing, among other keys, `'lambda'` (the estimated security level, in bits) and `'broken'` (a boolean indicating whether the scheme is considered broken at those parameters). For the example parameters above, the expected result is approximately `lambda ≈ 33.08` bits.



# Experiments

## Claim #1 — Security estimates for AGCD parameters from the literature

**Claim:** the estimator reproduces, for sets of `(gamma, eta, rho)` parameters taken from AGCD schemes proposed in the literature (DGHV, the Cheon et al. homomorphic scheme, Benarroch–Lepoint), the security level in bits reported/discussed in the paper.

**How to run:**

1. Inside the `sage` shell, at the repository root, start Jupyter from SageMath:

   ```bash
   sage -n jupyterlab
   ```

2. Open `agcd_test_notebook.ipynb`.

3. Run the cells in the `## Hilder params` section (parameters `gamma=680, eta=105, rho=100`), `## Cheon revisiting params` (parameters `gamma=157, eta=128, rho=121` and `gamma=1024, eta=128, rho=93`), and `## Benarroch and Lepoint` (parameters `gamma=175, eta=150, rho=110`), in that order.

   Alternatively, without Jupyter, the same results can be obtained directly in the `sage` shell:

   ```python
   from estimator import AGCD
   from estimator.reduction import RC

   AGCD.estimate(AGCD.Parameters(gamma=680, eta=105, rho=100))
   AGCD.estimate(AGCD.Parameters(gamma=157, eta=128, rho=121), red_cost_model=RC.MATZOV)
   AGCD.estimate(AGCD.Parameters(gamma=1024, eta=128, rho=93), red_cost_model=RC.ELHL26)
   AGCD.estimate(AGCD.Parameters(gamma=175, eta=150, rho=110))
   ```

**Relevant flags/configuration:** the `red_cost_model` parameter (e.g., `RC.MATZOV`, `RC.ELHL26`) selects the lattice-reduction cost model used in the estimate; use the same model cited in the paper for each parameter set.

**Expected time:** a few seconds per parameter set (all four together: under 1 minute).

**Expected resources:** under 500 MB of RAM; no significant disk usage.

**Expected result:** for each parameter set, the dictionary returned by `AGCD.estimate(...)` contains, under the `'final'` key, the `'lambda'` value (security level in bits) matching what is presented in Table/Section [X] of the paper. [FILL IN: exact reference numerical values.]

## Claim #2 — Empirical validation of the orthogonal lattice attack

**Claim:** the orthogonal lattice attack implemented in `ol_attack.sage` successfully recovers the secret divisor `p` from synthetic AGCD samples, on reduced-size instances, consistent with the feasibility predicted by the theoretical estimator in Claim #1.

**How to run:**

1. Inside the `sage` shell, at the repository root, run:

   ```bash
   sage ol_attack.sage
   ```

   This runs the `test()` function defined in the script itself, which by default uses `gamma=157, eta=128, rho=121, n=25` samples and `block_size=25` for BKZ.

2. To test other parameters, edit the `test()` call at the end of the file, or import the functions manually in the `sage` shell:

   ```python
   load("ol_attack.sage")
   p, samples = generate_agcd_samples(gamma=157, eta=128, rho=121, n=25)
   recover_p(samples, rho=121, block_size=25)
   ```

**Relevant flags/configuration:** `block_size` controls the BKZ blocksize used in lattice reduction; `n` controls the number of AGCD samples generated (should be compatible with the optimal `n` estimated in Claim #1 for the same `gamma, eta, rho`).

**Expected time:** approximately 10–60 seconds for the default parameters (`n=25`, `block_size=25`) on commodity hardware; can increase significantly for larger blocksizes.

**Expected resources:** under 1 GB of RAM for the default parameters.

**Expected result:** the script's output prints `[SUCCESS] Recovered a divisor of p (likely p itself).`, confirming that the candidate returned by the attack divides the generated secret `p` — empirically validating the attack's feasibility for the chosen parameters, as discussed in the paper. [FILL IN: reference to the specific Section/Table of the paper, if applicable.]

# LICENSE

This project is a fork of the [Lattice Estimator](https://github.com/malb/lattice-estimator), distributed under the [**LGPLv3+ (GNU Lesser General Public License v3.0 or later)**](https://www.gnu.org/licenses/lgpl-3.0.en.html), inherited from the original project.
