# ==============================
# Define ODE model & best fit parameters
# ==============================

bloodModel = "Combined HSC & MPP Self-Renewal"; 
include(joinpath(@__DIR__, "..", "Models", bloodModel * ".jl"));

df = CSV.read(joinpath(@__DIR__, "..", "Results", bloodModel, "WT ASA.csv"), DataFrame);

p = Matrix(df[:,4:end]); CF = df[:,3];
pbest = GetBestParam(p,CF)[1];

# ==============================
# Lower and upper bounds generation from data
# ==============================

t = time_d; IC = IC_WT;
num_points = 169;

#Individual level data
t_data = Float64.(unique(t));
idxs = (i28,i56,i84,i112);
cell_data = cellCount_WT;

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

function solve_ode(theta, time, IC, num_points)

    t_start = minimum(time); t_end = maximum(time)
    tspan = (t_start, t_end)
    stime = range(t_start, t_end, length=num_points)

    prob = ODEProblem(bloodODE!, IC, tspan, theta)
    sol = solve(prob, Tsit5(), saveat = stime, reltol=1e-6, abstol=1e-6);
    
    return sol
end

solbest = solve_ode(pbest, t, IC, num_points);
t_sim = solbest.t;

# ==============================
# Trajectories
# ==============================

df_p = CSV.read(joinpath(@__DIR__,"..", "Results", "Accepted Parameters (WT).csv"), DataFrame);

pop_names = namePops;
c_WT = ["#002856","#2B5A94","#97A9BD"];
title_WT = "H3-WT";

p = Matrix(df_p[:,3:end]);
acc_traj = size(p,1);

sols = Any[];
push!(sols, solbest);
for k in 2:acc_traj
    sol = solve_ode(p[k,:], t, IC, num_points);
    push!(sols, sol)
end
t_sim = sols[1].t;

# ==============================
# Figure S2 - WT
# ==============================

# H3-WT
nPops = size(namePops,1);

for i in 1:nPops

    f_WT = plot(vcat(
        [scatter(x = t_sim, y = sols[k][i, :], mode = "lines", name = "X", line = attr(width = 6, color = c_WT[3]), showlegend = false) for k in eachindex(sols)],
        [scatter(x=t_sim, y=solbest[i,:], mode="lines", name="X", line=attr(width=6, color=c_WT[1]), showlegend=false),
        scatter(x=t_data, y=lower_bounds[i,:], mode="markers", marker=attr(size=15, color="black"), showlegend=false),
        scatter(x=t_data, y=upper_bounds[i,:], mode="markers", marker=attr(size=15, color="black"), showlegend=false)]
    ))

    if namePops_WT[i] == "HSC1 + HSC2"
        width_f=750; height_f=375;
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
        width_f=375; height_f=375;

        if namePops_WT[i] == "MEP"
            xaxis_attr = attr(
                showticklabels = false,
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
                showticklabels = false,
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
    end

    f_WT.plot.layout.xaxis = xaxis_attr;
    f_WT.plot.layout.yaxis = yaxis_attr;

    relayout!(f_WT, 
        showlegend = false, 
        template = craig_lab_template,
        title=nothing,
        margin= attr(l=55, r=15, t=45, b=45), 
        width=width_f, height=height_f,
        annotations = [
            attr(
                text = namePops_WT[i],
                x = 0.5,
                y = 1.02,
                xref = "paper",
                yref = "paper",
                xanchor = "center",
                yanchor = "bottom",
                showarrow = false,
                font = attr(size = 32)
            )
        ]
    );

    display(f_WT)

end



