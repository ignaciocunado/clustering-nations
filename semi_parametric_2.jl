using Distributions
using Turing
using Turing.RandomMeasures
using StatsPlots
using DataFrames
using CSV
using LinearAlgebra

# Declare data path and import CSV into new DataFrame
path = "data/csv/OECD_data.csv";
df = DataFrame(CSV.File(path));
df.PRIVATESCH = map(x -> x == "private" ? 1 : 0, df.PRIVATESCH)


# Using Y_MATH1_rate as out Ys and the rest as covariates
# X = df[:, vcat(2:5, 7:9)]
X = df[:, 2:9]
Y = df[:, [2,11]]
grouped_X = groupby(X, :CNT) # Access school t in country i via grouped_X[i][t, :]
grouped_Y = groupby(Y, :CNT)

# Find some of the required data
n = length(Set(X[:, 2]))
n_i = combine(groupby(X, :CNT), nrow => :NumSchools)[:, 2]

# Declare the model
@model function hierarchical_dp_model(X, y, grouped_X, n, ni)
    # Hyperparameters
    M = 1.0  # DP concentration parameter (α)
    μ0 = 0.0
    σ0 = 10.0
    σ_β = 10.0
    
    # Dimensions
    p = size(X, 2) - 1
    
    # Prior for β
    β ~ MvNormal(zeros(p), σ_β^2 * I)
    
    # Create Dirichlet Process
    dp = DirichletProcess(M)
    G0 = Normal(μ0, σ0^2)
    
    # Initialize cluster assignments and effects
    z = tzeros(Int, n)  # Cluster assignments for each school
    c = tzeros(Float64, 0)  # Cluster-specific effects
    
    # Sample cluster assignments and effects using CRP
    for i in 1:n
        # Get current number of clusters (max value of z)
        K = maximum(z)
        nk = Vector{Int}(map(k -> sum(z .== k), 1:K))  # Count of schools in each cluster
        
        # Compute the probabilities of assigning to each cluster or a new cluster
        probs = Float64[]
        for k in 1:K
            push!(probs, nk[k] / (n - 1 + M))  # Probability of assigning to an existing cluster
        end
        
        # Probability of assigning to a new cluster
        push!(probs, M / (n - 1 + M))  # New cluster
        
        # Normalize the probabilities so that they sum to 1
        probs /= sum(probs)
        
        # Sample the cluster assignment z[i] based on computed probabilities
        z[i] ~ Categorical(probs)
        
        # If a new cluster is selected (i.e., z[i] > K), initialize a new cluster effect
        if z[i] > K
            push!(c, 0.0)  # New cluster effect
        end
        
        # Likelihood for schools in country i
        for t in 1:ni[i]
            λ = exp(dot(Vector(grouped_X[i][t, Not(:CNT)]), β) + c[z[i]] + 
                    log(grouped_X[i][t, :SCH_TESTED]))
            y[i][t, :Y_MATH1_rate] ~ Poisson(λ)
        end
    end
end


# Run sampling with initial parameters
chain = sample(
    hierarchical_dp_model(X, grouped_Y, grouped_X, n, n_i),
    SMC(),
    1000,
)

plot(chain)