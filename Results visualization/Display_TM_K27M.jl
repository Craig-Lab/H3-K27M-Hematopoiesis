
# ==============================
# Define ODE model & best fit parameters
# ==============================

bloodModel = "Transiting HSCs & MPP Self-Renewal"; 
include(bloodModel*".jl"); 

csvpath = pwd()*"/"*bloodModel*"/";
df = CSV.read(csvpath*"K27M ASA.csv", DataFrame);

p = Matrix(df[:,4:end]); CF = df[:,3];
pbest = GetBestParam(p,CF)[1];

fig_K27M = [time_d,cellCount_K27M,avgCell_K27M,stdCell_K27M];

# ==============================
# Lower and upper bounds generation from data
# ==============================

t = time_d; IC = IC_K27M;
num_points = 169;

#Individual level data
t_data = Float64.(unique(t));
idxs = (i28,i56,i84,i112);
cell_data = cellCount_K27M;

lower_bounds = zeros(Int64,nPops,4); 
upper_bounds = zeros(Int64,nPops,4); 

for i = 1:nPops
    # Extract datapoint at four the timepoints
    data_block  = @view cell_data[:,i]
    d1 = data_block[i28];
    d2 = data_block[i56];
    d3 = data_block[i84];
    d4 = data_block[i112];

    data_mat = hcat(d1, d2, d3, d4);
    lower_bounds[i,:] = minimum(data_mat, dims=1)
    upper_bounds[i,:] = maximum(data_mat, dims=1)
end

# ==============================
# ODE solver (used in the loop)
# ==============================

function solve_ode(theta, time, IC, num_points,i)

    t_start = minimum(time); t_end = maximum(time)
    tspan = (t_start, t_end)
    stime = range(t_start, t_end, length=num_points)

    prob = ODEProblem(bloodODE!, IC, tspan, theta)
    sol = solve(prob, Tsit5(), saveat = stime, reltol=1e-6, abstol=1e-6);
    
    return sol
end

sol = solve_ode(pbest, t, IC, num_points,1);
t_sim = sol.t;

# ==============================
# Trajectories
# ==============================

df_p = CSV.read("Accepted Parameters (K27M).csv", DataFrame);

pop_names = namePops;
c_K27M = ["#a0312b","#D44138","#DEBDBC"];
title_K27M = "H3-K27M";

p = Matrix(df_p[:,3:end]);
acc_traj = size(p,1);

sols = Any[];
for k in 2:acc_traj
    sol = solve_ode(p[k,:], t, IC, num_points,1);
    push!(sols, sol)
end
t_sim = sols[1].t;

solbest = solve_ode(pbest, t, IC, num_points,1);

pmean = vec(mean(p,dims = 1));
solmean = solve_ode(pmean, t, IC, num_points,1);

# ==============================
# 95% CI
# ==============================

# Find 95% confidence intervals
CI95_K27M = Vector{Matrix{Float64}}(undef,nPops);

for k in 1:nPops
    x = zeros(size(sols,1),num_points);
    CI95 = zeros(num_points,2);

    for j in eachindex(sols)
        x[j,:] = permutedims(sols[j][k,:])
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
    CI95_upper = Float64[]
);

for k in 1:nPops
    for i in 1:num_points

        push!(
            df_CI95_K27M,
            (
                pop_names[k],
                t_sim[i],
                CI95_K27M[k][i, 1],
                CI95_K27M[k][i, 2]
            )
        )

    end
end

# CSV.write("Solution CI95 K27M.csv", df_CI95_K27M);

# ==============================
# Figure S2 - K27M
# ==============================

# H3-K27M
nPops = size(namePops,1);

for i in 1:nPops

    f_K27M = plot(vcat(
        [scatter(x = t_sim, y = sols[k][i, :], mode = "lines", name = "X", line = attr(width = 6, color = c_K27M[3]), showlegend = false) for k in eachindex(sols)],
        [scatter(x=t_sim, y=solbest[i,:], mode="lines", name="X", line=attr(width=6, color=c_K27M[1]), showlegend=false),
        scatter(x=t_data, y=lower_bounds[i,:], mode="markers", marker=attr(size=15, color="black"), showlegend=false),
        scatter(x=t_data, y=upper_bounds[i,:], mode="markers", marker=attr(size=15, color="black"), showlegend=false)]
    ))

    width_f=375; height_f=375;

    if namePops_K27M[i] == "HSC1" || namePops_K27M[i] == "MEP"
        xaxis_attr = attr(
            tickvals = [28, 56, 84, 112],
            range = [20, 120],
            tickfont = attr(size = 28)
        );
        yaxis_attr = attr(
            type = "log",
            range = [2.7,7],
            tickmode = "array",
            tickvals = [exp10(3), exp10(5), exp10(7)],
            ticktext = ["10<sup>3</sup>", "10<sup>5</sup>", "10<sup>7</sup>"],
            tickfont = attr(size = 28)
        );
    else
        xaxis_attr = attr(
            tickvals = [28, 56, 84, 112],
            range = [20, 120],
            tickfont = attr(size = 28)
        );
        yaxis_attr = attr(
            type = "log",
            range = [2.7,7],
            tickmode = "array",
            tickvals = [exp10(3), exp10(5), exp10(7)],
            ticktext = ["10<sup>3</sup>", "10<sup>5</sup>", "10<sup>7</sup>"],
            tickfont = attr(size = 28),
            showticklabels = false
        );
    end

    f_K27M.plot.layout.xaxis = xaxis_attr;
    f_K27M.plot.layout.yaxis = yaxis_attr;

    annotations = if namePops_K27M[i] == "HSC1" || namePops_K27M[i] == "HSC2"
        [
            attr(
                text = namePops_K27M[i],
                x = 0.5,
                y = 1.0,
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = false,
                font = attr(size = 32)
            )
        ]
    else
        []
    end

    relayout!(f_K27M, 
        showlegend = false, 
        template = craig_lab_template,
        title=nothing,
        margin= attr(l=55, r=15, t=45, b=45), #attr(l=80, r=50, t=80, b=80), 
        width=width_f, height=height_f,
        annotations = annotations
    );

    display(f_K27M)
    figpath = pwd()*"/"*bloodModel*"/Figures/";
    # savefig(f_K27M,figpath*namePops_K27M[i]*" K27M - Trajectories.svg";width=width_f, height=height_f);

end