using Distributions
using Turing
using StatsPlots
using DataFrames
using CSV
using LinearAlgebra
using RCall


# Declare data path and import CSV into new DataFrame
path = "data/csv/OECD_data.csv"
df = DataFrame(CSV.File(path))
df.PRIVATESCH = map(x -> x == "private" ? 1 : 0, df.PRIVATESCH)


# Using Y_MATH1_rate as out Ys and the rest as covariates
X = df[:, vcat(2:6, 8:9)]
Y = df[:, [2,11]]
grouped_X = groupby(X, :CNT) # Access school t in country i via grouped_X[i][t, :]
grouped_Y = groupby(Y, :CNT)

# Find some of the required data
n = length(Set(X[:, 1]))
n_i = combine(groupby(X, :CNT), nrow => :NumSchools)[:, 2]


# Declare the model
@model function semi_parametric_approach(x, y, n_i, n, grouped_X, grouped_Y)

    # Hyperparameters
    sigma_b = 2.0
    mu_0 = 0
    sigma_0 = 2.0
    M = 5.0
    K = 20  # Truncation level for stick-breaking process

    β ~ MvNormal(zeros(size(x, 2) - 1), sigma_b ^ 2 * I)
    
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
        μ[k] ~ Normal(mu_0, sigma_0 ^ 2)
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


chain = sample(semi_parametric_approach(X, Y, n_i, n, grouped_X, grouped_Y), SMC(), 3000, discard=1000)

pyplot()
plot(chain)
savefig("semi_parametric_approach/chains.png")

chain_df = DataFrame(chain)

clus_d = chain_df[:, Cols(r"z_i")]
clus_d_int = mapcols(col -> Int.(col), clus_d)
clus = Matrix(clus_d_int)

@rput clus
R"library(salso)"
R"bestclust = salso(clus, loss=binder())"
@rget bestclust

CSV.write("data/clusters/clusters.csv", DataFrame(country = group_names = [key.CNT for key in keys(grouped_X)], cluster = bestclust))