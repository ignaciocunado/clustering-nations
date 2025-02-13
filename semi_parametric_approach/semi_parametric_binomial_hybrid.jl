using Distributions
using Turing
using StatsPlots
using DataFrames
using CSV
using LinearAlgebra
using RCall
using PyPlot
using Conda
using StatsFuns
using Statistics

#dataset
df = CSV.read("data/csv/OECD_data.csv", DataFrame)
df.PRIVATESCH = map(x -> x == "private" ? 1 : 0, df.PRIVATESCH)


# Using Y_MATH1_rate as out Ys and the rest as covariates
X = df[:, vcat(2:6, 9)]
numeric_cols = setdiff(names(X), ["CNT", "SCH_TESTED", "mean_ESCS_std"])
for col in numeric_cols
    X[!, col] = Float64.(X[!, col])  # Use `!` to modify in place
end

for col in numeric_cols
    X[!, col] = (X[!, col] .- mean(X[!, col])) ./ std(X[!, col])
end

Y = df[:, [2,10]]
grouped_X = groupby(X, :CNT) # Access school t in country i via grouped_X[i][t, :]
grouped_Y = groupby(Y, :CNT)

# Find some of the required data
n = length(Set(X[:, 1]))
n_i = combine(groupby(X, :CNT), nrow => :NumSchools)[:, 2]


# Declare the model
@model function semi_parametric_approach(x, y, n_i, n, grouped_X, grouped_Y)

    # Hyperparameters
    sigma_b = 5.0
    mu_0 = 0
    sigma_0 = 5.0
    M = 2.0
    K = 42  # Truncation level for stick-breaking process

    β ~ MvNormal(zeros(size(x, 2) - 2), sigma_b ^ 2 * I)
    
    # Stick-breaking process
    v = Vector{Real}(undef, K)
    for k in 1:K
        v[k] ~ Beta(1, M)
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
            trials = grouped_X[i][t, :SCH_TESTED]
            p = logistic(dot(Array(grouped_X[i][t, Not(:CNT, :SCH_TESTED)]), β) + μ[z_i[i]])
            grouped_Y[i][t, :Y_MATH1] ~ Binomial(trials, p)
        end
    end
end



sampler_gibbs = Gibbs(
    PG(15, :z_i),  # Particle Gibbs for the entire z_i array
    MH(:v),        # MH for the entire v vector
    MH(:μ),         # MH for the entire μ vector
    NUTS()
)

# Build the model
model = semi_parametric_approach(X, Y, n_i, n, grouped_X, grouped_Y)

# Run the sampling
chain = sample(
    model,
    sampler_gibbs,
    1000,          # total number of samples
    discard= 500   # burn-in
)

pyplot()
Plots.plot(chain)
Plots.savefig("chains_binomial_mix.png")

chain_df = DataFrame(chain)

clus_d = chain_df[:, Cols(r"z_i")]
clus_d_int = mapcols(col -> Int.(col), clus_d)
clus = Matrix(clus_d_int)

@rput clus
R"library(salso)"
R"bestclust = salso(clus, loss=VI())"
@rget bestclust

CSV.write("data/clusters/clusters_binomial_mix.csv", DataFrame(country = group_names = [key.CNT for key in keys(grouped_X)], cluster = bestclust))


# @model function retrieve_coefficients(x, y, n_i, n, grouped_X, grouped_Y, z_i)

#     # Hyperparameters
#     sigma_b = 6.0
#     mu_0 = 0
#     sigma_0 = 6.0
#     M = 1.0
#     K = length(z_i)  # Truncation level for stick-breaking process

#     β ~ MvNormal(zeros(size(x, 2) - 2), sigma_b ^ 2 * I)
    
#     # Stick-breaking process
#     v = Vector{Float64}(undef, K)
#     for k in 1:K
#         v[k] ~ Beta(1, M)
#     end
    
#     # Calculate stick-breaking weights
#     w = Vector{Float64}(undef, K)
#     w[1] = v[1]
#     for k in 2:K
#         w[k] = v[k] * prod(1 .- v[1:k-1])
#     end

#     w /= sum(w)
    
#     # Cluster parameters
#     μ = Vector{Float64}(undef, K)
#     for k in 1:K
#         μ[k] ~ Normal(mu_0, sigma_0)
#     end

#     # Likelihood
#     for i in 1:n
#         # Sample cluster assignment for each country
#         for t in 1:n_i[i]
#             trials = grouped_X[i][t, :SCH_TESTED]
#             p = logistic(dot(Array(grouped_X[i][t, Not(:CNT, :SCH_TESTED)]), β) + μ[z_i[i]])
#             grouped_Y[i][t, :Y_MATH1] ~ Binomial(trials, p)
#         end
#     end
# end

# sampler_gibbs2 = Gibbs(
#     MH(:β),        # MH for the entire β vector
#     MH(:v),        # MH for the entire v vector
#     MH(:μ)         # MH for the entire μ vector
# )

# # Now sample the model using the bestclust from salso
# coefficients_chain = sample(
#     retrieve_coefficients(X, Y, n_i, n, grouped_X, grouped_Y, bestclust),
#     sampler_gibbs2,
#     100,
#     discard= 10
# )

# coefficients_chain_df = DataFrame(coefficients_chain)
# coefficcents = coefficients_chain_df[:, r"^μ"]