# Clustering Nations by Educational Performance

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
