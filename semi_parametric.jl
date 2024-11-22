using Turing;
using StatsPlots;
using Distributions;
using DataFrames;
using CSV;

# Declare data path and import CSV into new DataFrame
path = "data/csv/OECD_data.csv";
df = DataFrame(CSV.File(path));


# Using Y_MATH1_rate as out Ys and the rest as covariates
X = df[:, 1:9]
Y = df[:, [2,11]]
grouped_X = groupby(X, :CNT) # Access school t in country i via grouped_X[i][t, :]
grouped_Y = groupby(Y, :CNT)

# Find some of the required data
n = length(Set(X[:, 2]))
n_i = combine(groupby(X, :CNT), nrow => :NumSchools)[:, 2]

# Declare the model
@model function semi_parametric_approach(x,y)

    

end

chain = sample(semi_parametric_approach(X, Y), NUTS(), 10000)

plot(chain)