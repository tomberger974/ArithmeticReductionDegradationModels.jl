using Pkg
Pkg.activate(".")

using Revise
using ArithmeticReductionDegradationModels
import ArithmeticReductionDegradationModels as ARD

using DataFrames


# Create an empty DataFrame with the specified columns and types
maintenances = DataFrame(
    DATE = Float64[],
    TYPE = String[]
)

# Optionally, you can add rows to the DataFrame
push!(maintenances, (2.46, "P"))  # Example: Adding a row with sample data

# Create an empty DataFrame with the specified columns and types
degradations = DataFrame(
    DATE = Float64[],
    NB_MAINTENANCES = Int64[],
    VALUE = Float64[],
    TYPE = String[]
)

for row in [(.53, 0, .6, "BETWEEN"), (1.22, 0, 1.65, "BETWEEN"), (2.70, 1, 1.18, "BETWEEN"), (3.47, 1, 2.03, "BETWEEN" )]
   push!(degradations, row)  # Example: Adding a row with sample data
end

for row in [(.5, 0, .6, "BETWEEN"),(2.7, 1, 1.8, "BETWEEN")]
   push!(degradations, row)  # Example: Adding a row with sample data
end

degradationdata = DegradationData(maintenances, degradations)


μ1 = 1.
μ2 = 2.
σ1 = 1.
ρ1 = [.5]

wienerARD1_1 = WienerARD1(["P"], μ1, σ1, ρ1)
wienerARD∞_1 = WienerARD∞(["P"], μ1, σ1, ρ1)

degradationdata
fit_mle(wienerARD1_1, degradationdata)
fit_mle(wienerARD∞_1, degradationdata)

ARD.loglikelihood(wienerARD1_1, degradationdata)
loglikelihood(wienerARD∞_1, [degradationdata], [μ1, σ1, ρ1[1]])