# Clustering Nations by Educational Performance

## The Project

This project aims to cluster nations based on their educational performance.

Data is taken from the OECD Program for International Students Assessment. It includes questionnaires from schools and students.

Two main approaches have been considered:

### 1. Bayesian Semiparametric Approach

$$
\begin{aligned}
    y_{it} \mid \beta, b_i &\stackrel{\text{ind}}{\sim} \mathcal{P}(\exp(X_{it}^T\beta + b_i + \log(T_{it}))), \\
    \beta &\sim \mathcal{N}_p(0, \sigma_\beta^2 I_p), \\
    b_i \mid G &\stackrel{\text{iid}}{\sim} G, \\
    G &\sim DP(M, G_0), \\
    G_0 &\sim \mathcal{N}(\mu_0, \sigma_0^2).
\end{aligned}
$$


Where:
- $M$: Precision parameter, controlling the variability of the Dirichlet process
- $T_{it}$: Number of students in country *i* and school *t*
- $log(T_{it})$: Offset term to normalize the count data
- $y_{it}$: Number of low-achieving students
- $b_i$: Clustering component from the Dirichlet process, shared by subjects in the same cluster

### 2. Nested Dirichlet Process

## Julia

### Installing Julia on macOS

If you have Homebrew installed, run
```bash
brew install --cask julia
```
Otherwise, run
```bash
curl -fsSL https://install.julialang.org | sh
```

### Installing Julia on Windows

Install the latest version from the Microsoft Store by running
```bash
winget install julia -s msstore
```

### Installing Julia on Linux

Run
```bash
curl -fsSL https://install.julialang.org | sh
```

### Activating the environment

Once Julia is installed, open a terminal, navigate to the project root and run the `julia` command.

Now run the following commands

```julia
using Pkg
Pkg.activate("julia_environment")
Pkg.instantiate()
```

Exit via `exit()`.

### Running Julia Scripts

On your terminal, navigate to the project root and run

```bash
julia <fileName>.jl
```

Otherwise, you can run the files using VSCode or an IDE of your choice.
