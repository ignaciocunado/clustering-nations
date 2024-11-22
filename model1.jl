using Pkg
Pkg.add(["DataFrames", "CSV", "Distributions", "Turing", "StatsPlots"])
using CSV, DataFrames, Turing, Distributions

data = CSV.read("OECD_data.csv", DataFrame)
data.mean_ESCS_std = (data.mean_ESCS .- mean(data.mean_ESCS)) ./ std(data.mean_ESCS)


@model function poisson_model(y, X, offsets)
    β ~ MvNormal(zeros(size(X, 2)), 10.0)  # Prior for coefficients
    c ~ filldist(Normal(0, 1), size(y, 1))  # Random effects
    
    for i in 1:length(y)
        λ = exp(X[i, :]' * β + c[i] + offsets[i])
        y[i] ~ Poisson(λ)
    end
end


## we dont need this part of sampling (?):
 
y = data.sum_MATH1below
X = hcat(data.mean_ESCS_std, data.STRATIO)  # Example covariates
offsets = log.(data.SCH_TESTED)

model = poisson_model(y, X, offsets)
chain = sample(model, NUTS(), 1000)  # 1000 samples with NUTS sampler

using StatsPlots
plot(chain)

view(data)
