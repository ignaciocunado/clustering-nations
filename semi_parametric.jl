using Distributions
using Turing
using StatsPlots
using DataFrames
using CSV
using LinearAlgebra
# using Turing.RandomMeasures;

# Declare data path and import CSV into new DataFrame
path = "data/csv/OECD_data.csv"
df = DataFrame(CSV.File(path))
df.PRIVATESCH = map(x -> x == "private" ? 1 : 0, df.PRIVATESCH)


# Using Y_MATH1_rate as out Ys and the rest as covariates
X = df[:, 2:9]
Y = df[:, [2,11]]
grouped_X = groupby(X, :CNT) # Access school t in country i via grouped_X[i][t, :]
grouped_Y = groupby(Y, :CNT)

# Find some of the required data
n = length(Set(X[:, 2]))
n_i = combine(groupby(X, :CNT), nrow => :NumSchools)[:, 2]


# Declare the model
@model function semi_parametric_approach(x, y, n_i, grouped_X, grouped_Y)

    # Hyperparameters
    sigma_b = 10.0
    mu_0 = 0
    sigma_0 = 10.0
    M = 2.0
    alpha_d = 1.0  # Concentration parameter for DP
    K = 20  # Truncation level for stick-breaking process

    β ~ MvNormal(zeros(size(x, 2) - 1), sigma_b ^ 2 * I)
    
    # Stick-breaking process
    v = Vector{Real}(undef, K)
    for k in 1:K
        v[k] ~ Beta()
    end
    
    # Calculate stick-breaking weights
    w = Vector{Real}(undef, K)
    w[1] = v[1]
    for k in 2:K
        w[k] = v[k] * prod(1 .- v[1:k-1])
    end

    w /= sum(w)
    
    # Cluster parameters
    μ = Vector{Real}(undef, K)
    for k in 1:K
        μ[k] ~ Normal(mu_0, sigma_0)
    end

    # Declare z_i as a vector of categorical variables (one per country)
    z_i = Vector{Categorical}(undef, n)  # For each country, assign a cluster

    # Likelihood
    for i in 1:n
        # Sample cluster assignment for each country
        z_i[i] ~ Categorical(w)
        
        for t in 1:n_i[i]
            λ = exp(dot(Array(grouped_X[i][t, Not(:CNT)]), β) + μ[z_i[i]] + log(grouped_X[i][t, :SCH_TESTED]))
            grouped_Y[i][t, :Y_MATH1_rate] ~ Poisson(λ)
        end
    end
end


chain = sample(semi_parametric_approach(X, Y, n_i, grouped_X, grouped_Y), SMC(), 1000)

plot(chain)