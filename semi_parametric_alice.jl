using Distributions
using Turing
using StatsPlots
using DataFrames
using CSV
using LinearAlgebra
using Turing.RandomMeasures
using MCMCChains


# Data import and preprocessing (unchanged)
df = CSV.read("data/csv/OECD_data.csv", DataFrame)
df.PRIVATESCH = ifelse.(df.PRIVATESCH .== "private", 1, 0)
df = dropmissing(df)  # Remove rows with missing values

X = df[:, vcat(2:5, 7:9)]
Y = df[:, [2,11]]
S = df[:, vcat(2:6, 7:9)]

grouped_X = groupby(X, :CNT)
grouped_Y = groupby(Y, :CNT)
grouped_X = [select(group, Not(:CNT)) for group in grouped_X]
grouped_S = groupby(S, :CNT)
grouped_S = [select(group, Not(:CNT)) for group in grouped_S]

n = length(Set(X[:, 2]))
n_i = combine(groupby(X, :CNT), nrow => :NumSchools)[:, 2]

# Hyperparameters
σβ = 10.0
μ0 = 0.0
σ0 = 10.0
α = 1.0  # Concentration parameter for DP
N = length(grouped_X)
p = size(grouped_X[1], 2)

@model function dirichlet_process_model(X, y, T_it, α, μ0, σ0, σβ)
    # Prior for β
    β ~ MvNormal(zeros(p), σβ^2 * I(p))

    # Dirichlet Process using Chinese Restaurant Process
    rpm = DirichletProcess(α)
    
    # Cluster parameters
    μ = tzeros(Float64, 0)
    
    # Latent assignments
    z = tzeros(Int, N)
    
    for i in 1:N
        # Number of existing clusters
        K = maximum(z)
        nk = [sum(z .== k) for k in 1:K]
        
        # Draw latent assignment
        z[i] ~ ChineseRestaurantProcess(rpm, convert(Vector{Int64}, nk))
        
        # Create new cluster if necessary
        if z[i] > K
            push!(μ, 0.0)
            μ[z[i]] ~ Normal(μ0, σ0)
        end
        
        # Likelihood of the model
        for t in 1:n_i[i]
            covariate_vector = Array(X[i][t,:])
            λ = exp(dot(covariate_vector, β) + μ[z[i]] + log(T_it[i][t, :SCH_TESTED]))
            y[i][t, :Y_MATH1_rate] ~ Poisson(λ)
        end
    end
end

# Sampling
chain = sample(dirichlet_process_model(grouped_X, grouped_Y, grouped_S, α, μ0, σ0, σβ), SMC(), 100)

# Plot results
plot(chain)
