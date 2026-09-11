using CSV, DataFrames, Distributions, DifferentialEquations, LinearAlgebra, Random, Base.Threads

# Add the name of any of the following model for "bloodModel": 
# "Combined HSC", "Combined HSC & MPP Self-Renewal", "Combined HSC & MEP Bypass", "Combined HSC & MPP Self-Renewal & MEP Bypass"
# "Transiting HSCs", "Transiting HSCs & MPP Self-Renewal", "Transiting HSCs & MEP Bypass", "Transiting HSCs & MPP Self-Renewal & MEP Bypass"

bloodModel = "Combined HSC";

# Load data and model
include(joinpath(@__DIR__, "..", "Models", bloodModel * ".jl"));
include("ASAFunctions.jl");

# Build CostData objects
save_times_WT = Float64.(unique(time_WT));
idxs_WT = (i28_WT,i56_WT,i84_WT,i112_WT);
costdata_WT = CostData(cellCount_WT, save_times_WT, idxs_WT);
p0_WT = L_WT .+ (U_WT .- L_WT).*rand();

save_times_K27M = Float64.(unique(time_K27M));
idxs_K27M = (i28_K27M,i56_K27M,i84_K27M,i112_K27M);
costdata_K27M = CostData(cellCount_K27M, save_times_K27M, idxs_K27M);
p0_K27M = L_K27M .+ (U_K27M .- L_K27M).*rand();

#Create ODE problems
tspan_WT = (minimum(time_WT), maximum(time_WT));
prob_WT  = ODEProblem(bloodODE!, IC_WT, tspan_WT, p0_WT);

tspan_K27M = (minimum(time_K27M), maximum(time_K27M));
prob_K27M  = ODEProblem(bloodODE!, IC_K27M, tspan_K27M, p0_K27M);

struct ASAResult
    run::Int
    condition::Symbol
    CF::Float64
    p::Vector{Float64}
end

# Initialization of ASA
ASA_runs = 200;

jobs = vcat(
    [(r, :WT)   for r in 1:ASA_runs],
    [(r, :K27M) for r in 1:ASA_runs] 
);

results = Vector{ASAResult}(undef,length(jobs));

Threads.@threads :dynamic for i in eachindex(jobs)

    run_id, condition = jobs[i];
    seed = hash((run_id, condition));
    rng = MersenneTwister(seed);

    if condition == :WT
        IC       = IC_WT;
        costdata = costdata_WT;
        L, U     = L_WT, U_WT;
        tspan    = tspan_WT;
    elseif condition == :K27M
        IC       = IC_K27M;
        costdata = costdata_K27M;
        L, U     = L_K27M, U_K27M;
        tspan    = tspan_K27M;
    end

    p0 = L .+ (U .- L).*rand(rng,length(L));
    prob = ODEProblem(bloodODE!, IC, tspan, p0);
    sol = ASA(p0, L, U, prob, costdata; rng=rng);

    p_best = sol[1]; CF_best = sol[2];

    results[i] = ASAResult(run_id, condition, CF_best, p_best);

end

# Extract fields
runs       = getfield.(results, :run);
conditions = getfield.(results, :condition);
CFs        = getfield.(results, :CF);
P          = getfield.(results, :p); 

# Build DataFrame
Pmat = reduce(vcat, transpose.(P));
param_names = Symbol.(labs);

df = DataFrame(
    run       = runs,
    condition = conditions,
    CF        = CFs
);
for j in 1:N_param
    df[!, param_names[j]] = Pmat[:, j]
end

# Save WT and K27M DataFrames
df_WT   = filter(:condition => ==(:WT), df);
df_K27M = filter(:condition => ==(:K27M), df);

# CSV.write(joinpath(@__DIR__, "..", "Results", bloodModel, "WT ASA.csv"), df_WT);
# CSV.write(joinpath(@__DIR__, "..", "Results", bloodModel, "K27M ASA.csv"), df_K27M);