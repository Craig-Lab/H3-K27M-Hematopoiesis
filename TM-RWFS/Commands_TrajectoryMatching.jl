using Statistics, Distributions, DifferentialEquations, StatsBase, PosteriorStats
using CSV, DataFrames

function sec2time(dt)
    #Transforms time in seconds to hours:minutes:seconds
    time = mod(dt,3600);
    hrs  = Int((dt - time)/3600);
    mins = Int((time-mod(time,60))/60);
    secs = Int(round(mod(time,60)));
    println("Elapsed time (H:M:S): ",hrs,":",mins,":",secs);
end

function GetBestParam(p,CF)

    CF_best = minimum(CF);
    idx = getindex.(indexin(CF_best,CF),1);
    p_best = p[idx,:];

    bestfit = [p_best, CF_best];
    return bestfit
end

function solve_ode(prob_template_WT, theta, stime)

    remake(prob_template_WT, p = theta)

    sol = solve(remake(prob_template_WT, p = theta),
                 Tsit5(),
                 saveat = stime,
                 reltol = 1e-6,
                 abstol = 1e-6)
    return sol
end

function bounds_generation(cell_data, nP)

    lbs = zeros(Int64,nP,4); 
    ubs = zeros(Int64,nP,4); 

    for i = 1:nP
        # Extract datapoint at four the timepoints
        data_block  = @view cell_data[:,i]
        d1 = data_block[i28];
        d2 = data_block[i56];
        d3 = data_block[i84];
        d4 = data_block[i112];

        data_mat = hcat(d1, d2, d3, d4);
        lbs[i,:] = minimum(data_mat, dims=1)
        ubs[i,:] = maximum(data_mat, dims=1)
    end
    return lbs, ubs
end

# =================================================================================

num_points = 169;
t_start = minimum(time_d); t_end = maximum(time_d);
stime = range(t_start, t_end, length=num_points);

# =================================================================================
# H3-WT

bloodModel_WT = "Combined HSC & MPP Self-Renewal"; 
include(joinpath(@__DIR__, "..", "Models", bloodModel_WT * ".jl"));

df_WT = CSV.read(joinpath(@__DIR__, "..", "Results", bloodModel_WT, "WT ASA.csv"), DataFrame);

p_WT = Matrix(df_WT[:,4:end]); CF_WT = df_WT[:,3];
pbest_WT = GetBestParam(p_WT,CF_WT)[1];

nPops_WT = nPops; labs_WT = labs_WT;
namePops_WT = namePops;
celln_WT = cellCount_WT; 

# ==============================
# Initialization
# ==============================

# ODE template
prob_template_WT = ODEProblem(bloodODE!, IC_WT, (t_start, t_end), pbest_WT);
solbest_WT = solve_ode(prob_template_WT, pbest_WT, stime);

# Lower and upper bounds generation from data
lbs_WT, ubs_WT = bounds_generation(celln_WT, nPops_WT);

# Standard deviation assumptions
pSD = 0.05*ones(size(pbest_WT,1));

M = 1; # Number of repeated runs (use M = 1: the model is simple and converges quickly)
theta_current = copy(pbest_WT);

U_WT[findall(x -> x == "p<sub>CMP→GMP</sub>", labs_WT)[1]] = 0.2;
U_WT[findall(x -> x == "p<sub>GMP→M</sub>", labs_WT)[1]] = 0.1;
L_WT[findall(x -> x == "p<sub>CMP→MEP</sub>", labs_WT)[1]] = 1.0;
L_WT[findall(x -> x == "p<sub>EP-→M</sub>", labs_WT)[1]] = 0.1;

df_TM_WT = DataFrame(nb=Int[], run=Int[]);

for lab in labs_WT
    df_TM_WT[!, Symbol(lab)] = Float64[]
end

target_nb = 1000;

# ==============================
# TM-RWFS algorithm
# ==============================

t0 = time();

nb = 0;

for m in 1:M
    println("Iteration $m")

    while nb < target_nb

        theta_proposal = similar(theta_current)

        for k in eachindex(theta_current)
            theta_proposal[k] = rand(Normal(theta_current[k], abs(pSD[k])))
        end

        sol = solve_ode(prob_template_WT, theta_proposal, stime)

        s = Int[]

        idx = [1, 56, 112, 169];

        for ii in eachindex(idx)
            if all(lbs_WT[j,ii] <= sol[j, idx[ii]] <= ubs_WT[j,ii] for j in 1:nPops)
                push!(s,ii)
            end
        end

        if length(s) == 4 &&
            all(theta_proposal[j] > L_WT[j] for j in eachindex(pbest_WT)) &&
            all(theta_proposal[j] < U_WT[j] for j in eachindex(pbest_WT))

            theta_current .= theta_proposal
            nb += 1
            # println(nb)

            push!(df_TM_WT, (;
                nb = nb,
                run = m,
                (Symbol(labs_WT[j]) => theta_current[j] for j in eachindex(labs_WT))...
            )) 
        end
    end
    nb >= target_nb && break
end

sec2time(time() - t0)

# ==============================
# 95% Credible intervals
# ==============================

pacc_WT = Matrix(df_TM_WT[:,3:end]);
acc_traj = size(pacc_WT,1);

sols_WT = Any[];
push!(sols_WT, solbest_WT);

for k in 2:acc_traj
    prob = ODEProblem(bloodODE!, IC_WT, (t_start, t_end), pacc_WT[k,:]);
    sol = solve_ode(prob, pacc_WT[k,:], stime);
    push!(sols_WT, sol)
end

CI95_WT = Vector{Matrix{Float64}}(undef,nPops_WT);

for k in 1:nPops_WT
    x = zeros(size(sols_WT,1),num_points);
    CI95 = zeros(num_points,2);

    for j in eachindex(sols_WT)
        x[j,:] = permutedims(sols_WT[j][k,:])
    end

    for i in 1:num_points
        HDI = hdi(x[:,i]; prob=0.95);
        max_HDI = maximum(HDI);
        min_HDI = minimum(HDI);
        CI95[i,:] = [min_HDI, max_HDI];
    end
    CI95_WT[k] = CI95;
end

df_CI95_WT = DataFrame(
    Population = String[],
    Time = Float64[],
    CI95_lower = Float64[],
    CI95_upper = Float64[],
);

for k in 1:nPops_WT
    for i in 1:num_points

        push!(
            df_CI95_WT,
            (
                namePops_WT[k],
                stime[i],
                CI95_WT[k][i, 1],
                CI95_WT[k][i, 2]
            )
        )

    end
end

# ==============================
# Save results
# ==============================

CSV.write(joinpath(@__DIR__, "..", "Results", "Accepted Parameters (WT).csv"), df_TM_WT);
CSV.write(joinpath(@__DIR__, "..", "Results", "Solution CI95 WT.csv"), df_CI95_WT);

# =================================================================================
# H3-K27M
bloodModel_K27M = "Transiting HSCs & MPP Self-Renewal"; 
include(joinpath(@__DIR__, "..", "Models", bloodModel_K27M * ".jl"));

df_K27M = CSV.read(joinpath(@__DIR__, "..", "Results", bloodModel_K27M, "WT ASA.csv"), DataFrame);

p_K27M = Matrix(df_K27M[:,4:end]); CF_K27M = df_K27M[:,3];
pbest_K27M = GetBestParam(p_K27M,CF_K27M)[1];

nPops_K27M = nPops; labs_K27M = labs;
namePops_K27M = namePops;
celln_K27M = cellCount_K27M; 

# ==============================
# Initialization
# ==============================

# ODE template
prob_template_K27M = ODEProblem(bloodODE!, IC_K27M, (t_start, t_end), pbest_K27M);

# Lower and upper bounds generation from data
lbs_K27M, ubs_K27M = bounds_generation(celln_K27M, nPops_K27M);

# Standard deviation assumptions
pSD = 0.05*ones(size(pbest_K27M,1));

M = 1; # Number of repeated runs (use M = 1: the model is simple and converges quickly)

theta_current = copy(pbest_K27M);

U_K27M[findall(x -> x == "p<sub>CMP→GMP</sub>", labs_K27M)[1]] = 0.1;
U_K27M[findall(x -> x == "p<sub>GMP→M</sub>", labs_K27M)[1]] = 0.01;
U_K27M[findall(x -> x == "d<sub>MEP</sub>", labs_K27M)[1]] = 0.2;
U_K27M[findall(x -> x == "p<sub>EP-→M</sub>", labs_K27M)[1]] = 3;
U_K27M[findall(x -> x == "r<sub>HSC1→HSC2</sub>", labs_K27M)[1]] = 0.1;
U_K27M[findall(x -> x == "r<sub>HSC2→HSC1</sub>", labs_K27M)[1]] = 0.1;

df_TM_K27M = DataFrame(nb=Int[], run=Int[]);

for lab in labs_K27M
    df_TM_K27M[!, Symbol(lab)] = Float64[]
end

target_nb = 1000;

# ==============================
# TM-RWFS algorithm
# ==============================

t0 = time();

nb = 0;

for m in 1:M
    println("Iteration $m")

    while nb < target_nb

        theta_proposal = similar(theta_current)

        for k in eachindex(theta_current)
            theta_proposal[k] = rand(Normal(theta_current[k], abs(pSD[k])))
        end

        sol = solve_ode(prob_template_K27M, theta_proposal, stime)

        s = Int[]
        idx = [1, 56, 112, 169];

        for ii in eachindex(idx)
            if ii == 2
                if all(lbs_K27M[j,ii] <= sol[j, idx[ii]] <= ubs_K27M[j,ii] for j in 1:nPops if j != 2)
                    push!(s,ii)
                end
            else
                if all(lbs_K27M[j,ii] <= sol[j, idx[ii]] <= ubs_K27M[j,ii] for j in 1:nPops)
                    push!(s,ii)
                end
            end
        end

        if length(s) == 4 &&
            all(theta_proposal[j] > L_K27M[j] for j in eachindex(pbest_K27M)) &&
            all(theta_proposal[j] < U_K27M[j] for j in eachindex(pbest_K27M))

            theta_current .= theta_proposal
            nb += 1

            push!(df_TM_K27M, (;
                nb = nb,
                run = m,
                (Symbol(labs_K27M[j]) => theta_current[j] for j in eachindex(labs_K27M))...
            )) 
        end
    end

    nb >= target_nb && break
end

sec2time(time() - t0)

# ==============================
# 95% Credible intervals
# ==============================

pacc_K27M = Matrix(df_TM_K27M[:,3:end]);
acc_traj = size(pacc_K27M,1);

sols_K27M = Any[];
push!(sols_K27M, solbest_K27M);

for k in 2:acc_traj
    prob = ODEProblem(bloodODE!, IC_K27M, (t_start, t_end), pacc_K27M[k,:]);
    sol = solve_ode(prob, pacc_K27M[k,:], stime);
    push!(sols_K27M, sol)
end

CI95_K27M = Vector{Matrix{Float64}}(undef,nPops_K27M);

for k in 1:nPops_K27M
    x = zeros(size(sols_K27M,1),num_points);
    CI95 = zeros(num_points,2);

    for j in eachindex(sols_K27M)
        x[j,:] = permutedims(sols_K27M[j][k,:])
    end

    for i in 1:num_points
        HDI = hdi(x[:,i]; prob=0.95);
        max_HDI = maximum(HDI);
        min_HDI = minimum(HDI);
        CI95[i,:] = [min_HDI, max_HDI];
    end
    CI95_K27M[k] = CI95;
end

df_CI95_K27M = DataFrame(
    Population = String[],
    Time = Float64[],
    CI95_lower = Float64[],
    CI95_upper = Float64[],
);

for k in 1:nPops_K27M
    for i in 1:num_points

        push!(
            df_CI95_WT,
            (
                namePops_K27M[k],
                stime[i],
                CI95_K27M[k][i, 1],
                CI95_K27M[k][i, 2]
            )
        )

    end
end

# ==============================
# Save results
# ==============================

CSV.write(joinpath(@__DIR__, "..", "Results", "Accepted Parameters (K27M).csv"), df_TM_K27M);
CSV.write(joinpath(@__DIR__, "..", "Results", "Solution CI95 K27M.csv"), df_CI95_K27M);