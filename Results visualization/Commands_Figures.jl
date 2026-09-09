using PlotlyJS, CSV, DataFrames, DifferentialEquations, Statistics, StatsBase, HypothesisTests, PosteriorStats

include("FiguresFunctions.jl");
include("craig_lab_template.jl");

bloodModel_WT = "Combined HSC & MPP Self-Renewal"; 
bloodModel_K27M = "Transiting HSCs & MPP Self-Renewal";

# ---------------------------------------------------------------------- %
# Result for best fits (H3-WT)  

include(bloodModel_WT*".jl");
csvpath = pwd()*"/"*bloodModel_WT*"/";
df_WT = CSV.read(csvpath*"WT ASA.csv", DataFrame);

p_WT = Matrix(df_WT[:,4:end]);
CF_WT = df_WT[:,3];
ASA_WT = [p_WT,CF_WT];
pbest_WT = GetBestParam(p_WT,CF_WT)[1];

labs_WT = labs; namePops_WT = namePops; namePops_WT[1] = "HSC1 + HSC2";
fig_WT = [time_d,cellCount_WT,avgCell_WT,stdCell_WT];
nPops_WT = nPops;

end_idx = findall(x -> x == 112, time_d);
timeWT_end = time_d[end_idx];
cellCountWT_end = cellCount_WT[end_idx,:];

# Solution for best fit results
sol_WT = SolveBloodODE(pbest_WT,time_d,IC_WT);

# Correlation matrices (within the same model)
corr_matrix_WT = cor(p_WT);

# ---------------------------------------------------------------------- %
# Result for best fits (H3-K27M)  

include(bloodModel_K27M*".jl");
csvpath = pwd()*"/"*bloodModel_K27M*"/";
df_K27M = CSV.read(csvpath*"K27M ASA.csv", DataFrame);

p_K27M = Matrix(df_K27M[:,4:end]);
CF_K27M = df_K27M[:,3];
ASA_K27M = [p_K27M,CF_K27M];
pbest_K27M = GetBestParam(p_K27M,CF_K27M)[1];

labs_K27M = labs; namePops_K27M = namePops;
fig_K27M = [time_d,cellCount_K27M,avgCell_K27M,stdCell_K27M];
nPops_K27M = nPops;

end_idx = findall(x -> x == 112, time_d);
timeK27M_end = time_d[end_idx];
cellCountK27M_end = cellCount_K27M[end_idx,:];
avgK27M_end = avgCell_K27M[end,:];
stdK27M_end = stdCell_K27M[end,:];

# Solution for best fit results
sol_K27M = SolveBloodODE(pbest_K27M,time_d,IC_K27M);

# Correlation matrices (within the same model)
corr_matrix_K27M = cor(p_K27M);

# ---------------------------------------------------------------------- %
# Setting up the info

meanTime = unique(time_d);

com = intersect(labs_K27M, labs_WT);
uni_WT = setdiff(labs_WT, labs_K27M);
uni_K27M = setdiff(labs_K27M, labs_WT);

idx_WT   = [findfirst(==(x), labs_WT)   for x in com];
idx_K27M = [findfirst(==(x), labs_K27M) for x in com];
idx_unique_WT   = [findfirst(==(x), labs_WT)   for x in uni_WT];
idx_unique_K27M = [findfirst(==(x), labs_K27M) for x in uni_K27M];

common_WT   = p_WT[:,idx_WT];
common_K27M = p_K27M[:,idx_K27M];
unique_WT   = p_WT[:,idx_unique_WT];
unique_K27M = p_K27M[:,idx_unique_K27M];

c_WT = ["#002856","#2B5A94","#97A9BD"];
c_K27M = ["#a0312b","#D44138","#DEBDBC"];
title_WT = "H3-WT";
title_K27M = "H3-K27M";

## =============================
# Figure 2
# ==============================

# H3-WT
df_CI_WT = CSV.read("Solution CI95 WT.csv", DataFrame);
groups_WT = groupby(df_CI_WT,:Population);

meanTime = unique(time_d);

for i in 1:nPops_WT
    time_CI_WT = Vector{Float64}(groups_WT[i][!, :Time]);
    CI_low_WT = Vector{Float64}(groups_WT[i][!, :CI95_lower]);
    CI_up_WT = Vector{Float64}(groups_WT[i][!, :CI95_upper]);

    f_WT = plot([
        scatter(x = time_CI_WT, y = CI_up_WT, mode = "lines", line = attr(color = "transparent"), line_width = 6.0, name = "97.5%", showlegend = false),
        scatter(x = sol_WT.t, y = CI_low_WT, mode = "lines", fill = "tonexty", fillcolor = hex2rgba(c_WT[2],0.1), line = attr(color = "transparent"), line_width = 6.0, name = "2.5%", showlegend = false),
        scatter(x = time_d, y = cellCount_WT[:,i], mode = "markers", marker = attr(size = 15, color = c_WT[3]), name = "Data points", visible = true),
        scatter(x = meanTime, y = avgCell_WT[:,i], mode = "markers", marker = attr(size = 15, color = c_WT[1], symbol = 1), name = "Mean values"),
        scatter(x = sol_WT.t, y = sol_WT[i,:], mode = "lines", marker = attr(color = c_WT[2]), line_width = 6.0, name = "Solution", visible = true)
    ]);

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
    # figpath = pwd()*"/"*bloodModel_WT*"/Figures/";
    # savefig(f_WT,figpath*namePops_WT[i]*" WT.svg";width=width_f, height=height_f);

end

# H3-K27M
df_CI_K27M = CSV.read("Solution CI95 K27M.csv", DataFrame);
groups_K27M = groupby(df_CI_K27M,:Population);

for i in 1:nPops_K27M
    time_CI_K27M = Vector{Float64}(groups_K27M[i][!, :Time]);
    CI_low_K27M = Vector{Float64}(groups_K27M[i][!, :CI95_lower]);
    CI_up_K27M = Vector{Float64}(groups_K27M[i][!, :CI95_upper]);

    f_K27M = plot([
        scatter(x = time_CI_K27M, y = CI_up_K27M, mode = "lines", line = attr(color = "transparent"), line_width = 6.0, name = "97.5%", showlegend = false),
        scatter(x = sol_K27M.t, y = CI_low_K27M, mode = "lines", fill = "tonexty", fillcolor = hex2rgba(c_K27M[2],0.2), line = attr(color = "transparent"), line_width = 6.0, name = "2.5%", showlegend = false),
        scatter(x = time_d, y = cellCount_K27M[:,i], mode = "markers", marker = attr(size = 15, color = c_K27M[3]), name = "Data points", visible = true),
        scatter(x = meanTime, y = avgCell_K27M[:,i], mode = "markers", marker = attr(size = 15, color = c_K27M[1], symbol = 1), name = "Mean values"),
        scatter(x = sol_K27M.t, y = sol_K27M[i,:], mode = "lines", marker = attr(color = c_K27M[2]), line_width = 6.0, name = "Solution", visible = true)
    ]);

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
        margin= attr(l=55, r=15, t=45, b=45),
        width=width_f, height=height_f,
        annotations = annotations
    );

    display(f_K27M)
    # figpath = pwd()*"/"*bloodModel_K27M*"/Figures/";
    # savefig(f_K27M,figpath*namePops_K27M[i]*" K27M.svg";width=width_f, height=height_f);

end

## =============================
# Figure 3A
# ==============================

pbestc_WT = pbest_WT[idx_WT];
pbestc_K27M = pbest_K27M[idx_K27M];

N = size(common_WT,2);

for i in eachindex(com)
    xmin1 = minimum(common_WT[:,i]);
    xmax1 = maximum(common_WT[:,i]);
    xmin2 = minimum(common_K27M[:,i]);
    xmax2 = maximum(common_K27M[:,i]);

    xmin = max(0.0, min(xmin1, xmin2));
    xmax = max(xmax1, xmax2);

    xmin_pad, xmax_pad, ticks, step = nice_ticks(xmin, xmax);

    # determine decimal formatting for x axis
    decimals = tick_decimals(step);
    tick_text = tick_label.(ticks, decimals);

    # histogram binning
    nbins = 15;
    bin_size = (xmax_pad - xmin_pad) / nbins;

    yvals1 = histogram_y_extent(common_WT[:,i], xmin_pad, xmax_pad, bin_size);
    ymax = maximum(yvals1);
    yvals2 = histogram_y_extent(common_K27M[:,i], xmin_pad, xmax_pad, bin_size);
    ymax = max(ymax, maximum(yvals2));

    # Histogram
    his = plot([
        histogram(
            x = common_WT[:,i],
            xbins = attr(
                start = xmin_pad,
                stop  = xmax_pad,
                size  = bin_size
            ),
            showlegend = false,
            marker = attr(
                color = c_WT[1],
                opacity = 0.30,
                line = attr(width=2.0, color=c_WT[1])
            )
        ),
        histogram(
            x = common_K27M[:,i],
            xbins = attr(
                start = xmin_pad,
                stop  = xmax_pad,
                size  = bin_size
            ),  
            showlegend = false,
            marker = attr(
                color = c_K27M[1],
                opacity = 0.30,
                line = attr(width=2.0, color=c_K27M[1])
            )
        )
    ]);

    y_min_pad, y_max_pad, y_ticks, y_step = nice_ticks(0.0, ymax);
    y_decimals = tick_decimals(y_step);
    y_tick_text = tick_label.(y_ticks, y_decimals);

    xaxis_attr = attr(
        range = [xmin_pad, xmax_pad],
        tickmode = "array",
        tickvals = ticks,
        ticktext = tick_text,
        tickformat = "",
        automargin = true,
        tickfont = attr(size = 28),
        showgrid = false
    );

    yaxis_attr = attr(
        range = [y_min_pad, y_max_pad],
        tickmode = "array",
        tickvals = y_ticks,
        ticktext = y_tick_text,
        tickformat = "",
        tickfont = attr(size = 28)
    );

    his.plot.layout.xaxis = xaxis_attr;
    his.plot.layout.yaxis = yaxis_attr;

    width_f=375; height_f=375;

    annotations =[
        attr(
            text = com[i],
            x = 0.5,
            y = 1.0,
            xref = "paper",
            yref = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = false,
            font = attr(size = 32)
        )
    ];

    relayout!(
        his,
        template = craig_lab_template,
        annotations = annotations,
        showlegend = false,
        margin = attr(l=55, r=15, t=45, b=45),
        width = width_f,
        height = height_f,
        barmode = "overlay",
        autosize = false,
        shapes = [
            vline(
                pbestc_WT[i],
                line = attr(
                    color = c_WT[2],
                    width = 4.0,
                    dash = "dash"
                )
            ),
            vline(
                pbestc_K27M[i],
                line = attr(
                    color = c_K27M[2],
                    width = 4.0,
                    dash = "dash"
                )
            )
        ]
    );

    display(his)
    # figpath = pwd()*"/";
    # savefig(his,figpath*"Param $i.svg";width=width_f, height=height_f);
end

## =============================
# Figure 3B
# ==============================

desired_order = [
    "k<sub>MPP</sub>", 
    "p<sub>MPP→MLP</sub>", 
    "p<sub>MPP→CMP</sub>", 
    "p<sub>CMP→GMP</sub>", 
    "p<sub>CMP→MEP</sub>", 
    "p<sub>MEP→EP+</sub>", 
    "p<sub>EP+→EP-</sub>", 
    "p<sub>MLP→M</sub>", 
    "p<sub>GMP→M</sub>", 
    "p<sub>EP-→M</sub>",
    "d<sub>CMP</sub>", 
    "d<sub>MEP</sub>", 
    "d<sub>EP+</sub>"
];

index = Dict(l => i for (i, l) in enumerate(desired_order));
rc = (pbestc_K27M .- pbestc_WT)./abs.(pbestc_WT);

# enforce order AFTER masking
perm = sortperm(com, by = x -> get(index, x, Inf));

com = com[perm];
rc = rc[perm];
pbestc_WT  = pbestc_WT[perm];
pbestc_K27M = pbestc_K27M[perm];

trace = scatter(
    x = com,
    y = rc,
    mode = "markers",
    marker = attr(color=c_K27M[2], size=20),
    showlegend=false
);

lines = [
    scatter(
        x = [com[i], com[i]],
        y = [0, rc[i]],
        mode = "lines",
        line = attr(color="black", width=4),
        showlegend=false
    )
    for i in eachindex(com)
    if !isnan(pbestc_WT[i]) && !isnan(pbestc_K27M[i])
];

annotations = [
    attr(
        x = 0.00,
        y = -0.40,
        xref = "paper",
        yref = "paper",
        text = "<b>Proliferation</b>",
        showarrow = false,
        xanchor = "left",
        font = attr(
            size = 24,
            color = "#222222"
        )
    ),
    attr(
        x = 0.43,
        y = -0.40,
        xref = "paper",
        yref = "paper",
        text = "<b>Differentiation / Outflow</b>",
        showarrow = false,
        xanchor = "center",
        font = attr(
            size = 24,
            color = "#222222"
        )
    ),
    attr(
        x = 0.91,
        y = -0.40,
        xref = "paper",
        yref = "paper",
        text = "<b>Death</b>",
        showarrow = false,
        xanchor = "center",
        font = attr(
            size = 24,
            color = "#222222"
        )
    )
];

separator_shapes = [
    attr(
        type = "line",
        x0 = 0.5,
        x1 = 0.5,
        y0 = 0,
        y1 = 1,
        xref = "x",
        yref = "paper",
        line = attr(
            color = "#DDDDDD",
            width = 2
        ),
        layer = "below"
    ),
    attr(
        type = "line",
        x0 = 9.5,
        x1 = 9.5,
        y0 = 0,
        y1 = 1,
        xref = "x",
        yref = "paper",
        line = attr(
            color = "#DDDDDD",
            width = 2
        ),
        layer = "below"
    )
];

width_f  = 4.5 * 375;
height_f = 1.5 * 375;

layout = Layout(
    template = craig_lab_template,
    xaxis = attr(
        categoryorder = "array",
        categoryarray = desired_order,
        tickangle = 45,
        automargin = true,
        tickfont = attr(
            size = 28,
            color = "#222222"
        ),
        ticklabelposition = "outside"
    ),
    yaxis = attr(
        title = "",
        nticks = 5,
        tickfont = attr(
            size = 28,
            color = "#222222"
        ),
        zeroline = true,
        zerolinecolor = "#555555",
        zerolinewidth = 1.5
    ),
    margin = attr(
        l = 55,
        r = 15,
        t = 45,
        b = 150
    ),
    height = Int(round(height_f)),
    width  = Int(round(width_f)),
    annotations = annotations,
    shapes = separator_shapes
);

hcomp = plot(vcat(lines, [trace]...), layout);
display(hcomp)
# savefig(hcomp,pwd()*"/Parameter Comparison.svg"; width = Int(round(width_f)), height = Int(round(height_f)));

## =============================
# Figure S2
# ==============================

include("Display_TM_WT.jl");
include("Display_TM_K27M.jl");

## =============================
# Figure S3
# ==============================

jitter = 0.12;

df_tTest = DataFrame(
    Population = String[],
    pvalue = Float64[],
    tstatistic = Float64[],
    df = Float64[]
);

for i in 1:nPops_WT

    if i == 1

        test = UnequalVarianceTTest(cellCountWT_end[:,i], cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1]);
        pval = pvalue(test); tstat = test.t; dof = test.df;

        push!(
            df_tTest,
            (
                namePops_WT[i],
                pval,
                tstat,
                dof
            )
        )

        x_WT   = 1 .+ (rand(length(cellCountWT_end[:,i])) .- 0.5) .* 2 .* jitter;
        x_K27M = 2 .+ (rand(length(cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1])) .- 0.5) .* 2 .* jitter;

        p_WT = scatter(
            x = x_WT,
            y = cellCountWT_end[:,i],
            mode = "markers", 
            marker = attr(size = 15, color = c_WT[3]),
            hoverinfo = "skip",
            showlegend = false
        );

        p_K27M = scatter(
            x = x_K27M,
            y = cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1],
            mode = "markers", 
            marker = attr(size = 15, color = c_K27M[3]),
            hoverinfo = "skip",
            showlegend = false
        );

        # Means
        m_WT = mean(cellCountWT_end[:,i]);
        m_K27M = mean(cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1]);

        # Standard error
        se_WT = std(cellCountWT_end[:,i])  / sqrt(length(cellCountWT_end[:,i]));
        se_K27M = std(cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1]) / sqrt(length(cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1]));

        minval = minimum(vcat(cellCountWT_end[:,i],cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1],m_WT - se_WT,m_K27M - se_K27M));
        maxval = maximum(vcat(cellCountWT_end[:,i],cellCountK27M_end[:,i] .+ cellCountK27M_end[:,i+1],m_WT + se_WT,m_K27M + se_K27M));
        
        ymin_pad, ymax_pad, y_ticks, y_step = nice_ticks(minval, maxval, n_ticks = 5);

    else

        test = UnequalVarianceTTest(cellCountWT_end[:,i], cellCountK27M_end[:,i+1]);
        pval = pvalue(test); tstat = test.t; dof = test.df;

        push!(
            df_tTest,
            (
                namePops_WT[i],
                pval,
                tstat,
                dof
            )
        )

        x_WT   = 1 .+ (rand(length(cellCountWT_end[:,i])) .- 0.5) .* 2 .* jitter;
        x_K27M = 2 .+ (rand(length(cellCountK27M_end[:,i+1])) .- 0.5) .* 2 .* jitter;

        p_WT = scatter(
            x = x_WT,
            y = cellCountWT_end[:,i],
            mode = "markers", 
            marker = attr(size = 15, color = c_WT[3]),
            hoverinfo = "skip",
            showlegend = false
        );

        p_K27M = scatter(
            x = x_K27M,
            y = cellCountK27M_end[:,i+1],
            mode = "markers", 
            marker = attr(size = 15, color = c_K27M[3]),
            hoverinfo = "skip",
            showlegend = false
        );

        # Means
        m_WT = mean(cellCountWT_end[:,i]);
        m_K27M = mean(cellCountK27M_end[:,i+1]);

        # Standard error
        se_WT = std(cellCountWT_end[:,i])  / sqrt(length(cellCountWT_end[:,i]));
        se_K27M = std(cellCountK27M_end[:,i+1]) / sqrt(length(cellCountK27M_end[:,i+1]));

        minval = minimum(vcat(cellCountWT_end[:,i],cellCountK27M_end[:,i+1],m_WT - se_WT,m_K27M - se_K27M));
        maxval = maximum(vcat(cellCountWT_end[:,i],cellCountK27M_end[:,i+1],m_WT + se_WT,m_K27M + se_K27M));
        
        ymin_pad, ymax_pad, y_ticks, y_step = nice_ticks(minval, maxval, n_ticks = 5);
    end

    power = floor(log10(y_ticks[end]));
    expmax = 10^power;

    # Summary scatters
    summary_WT = summary_trace(1, m_WT, se_WT);
    summary_K27M = summary_trace(2, m_K27M, se_K27M);

    exp_y = attr(
        text="×10<sup>$(Int.(power))</sup>",
        x=0, y=1.0,
        xref = "paper", yref = "paper",
        xanchor="left", yanchor="bottom",
        showarrow=false,
        font=attr(size=24)
    )

    layout_f = Layout(
        template = craig_lab_template,
        width = 375, height = 375,
        annotations = [exp_y,
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
            )],
        xaxis = attr(
            automargin = false,
            tickmode = "array",
            tickvals = [1, 2],
            ticktext = ["H3-WT", "H3-K27M"],
            range = [0.5, 2.5],
            showgrid = false,
            zeroline = false,
            ticks = "",
            tickfont = attr(size=28)
        ),
        yaxis = attr(
            automargin = false,
            title = "",
            tickvals = collect(y_ticks),
            ticktext = string.(round.(y_ticks ./ expmax ; digits = 1)),
            range = [ymin_pad, ymax_pad],
            showgrid = false,
            zeroline = false,
            tickfont = attr(size=28),
            titlefont = attr(size=32)
        ),
        showlegend = false,
        margin = attr(l=55, r=15, t=45, b=45)
    );

    common_WT = plot([p_WT, p_K27M, summary_WT..., summary_K27M...],layout_f);
    display(common_WT)
    # figpath = pwd()*"/";
    # savefig(common_WT,figpath*namePops_WT[i]*" - End time_d point.svg";width=375, height=375);
end

# CSV.write("tTest (unequal variance).csv", df_tTest);

## =============================
# Figure S4A
# ==============================

c = ["#000000","E25140"];
time_sol = sol_WT.t;

for i in 1:nPops_WT

    if i == 1

        norm = (sol_K27M[i,:] .+ sol_K27M[i+1,:])./sol_WT[i,:];

        f_norm = plot([
            scatter(x = time_sol, y = norm, mode = "lines", marker = attr(color = c[1]), line_width = 6.0, name = "Solution", visible = true),
            scatter(x = 21:0.5:119, y = ones(size(21:0.5:119)), mode = "lines", marker = attr(color = c[2]), line_width = 6.0, line_dash = "dot", name = "Solution", visible = true)
        ])
    else

        norm = sol_K27M[i+1,:]./sol_WT[i,:];

        f_norm = plot([
            scatter(x = time_sol, y = norm, mode = "lines", marker = attr(color = c[1]), line_width = 6.0, name = "Solution", visible = true),
            scatter(x = 21:0.5:119, y = ones(size(21:0.5:119)), mode = "lines", marker = attr(color = c[2]), line_width = 6.0, line_dash = "dot", name = "Solution", visible = true)
        ])
    end

    ymin = maximum(vcat(0.0,minimum(norm)));
    ymax = maximum(norm);
    ymin_pad, ymax_pad, y_ticks, y_step = nice_ticks(ymin, ymax, n_ticks = 5);
        
    width_f=375; height_f=375;

    index = findfirst(==(namePops_WT[i]), namePops_WT);

    if index <= 3
        xaxis_attr = attr(
            showticklabels = false,
            tickvals = [28, 56, 84, 112],
            range = [20, 120],
            tickfont = attr(size = 28)
        );
        yaxis_attr = attr(
            showticklabels = true,
            range = [ymin_pad, ymax_pad],
            tickmode = "array",
            tickvals = y_ticks,
            tickfont = attr(size = 28)
        );
    else
        xaxis_attr = attr(
            showticklabels = true,
            tickvals = [28, 56, 84, 112],
            range = [20, 120],
            tickfont = attr(size = 28)
        );
        yaxis_attr = attr(
            showticklabels = true,
            range = [ymin_pad, ymax_pad],
            tickmode = "array",
            tickvals = y_ticks,
            tickfont = attr(size = 28)
        );
    end

    f_norm.plot.layout.xaxis = xaxis_attr;
    f_norm.plot.layout.yaxis = yaxis_attr;

    relayout!(f_norm, 
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

    display(f_norm)
    # figpath = pwd()*"/";
    # savefig(f_norm,figpath*namePops_WT[i]*" - Normalized WT.svg";width=width_f, height=height_f);
end

## =============================
# Figure S4B
# ==============================
    
time_sol = sol_WT.t;

for i in 1:nPops_WT

    if i == 1

        norm_WT = sol_WT[i,:]./sol_WT[i,1];
        norm_K27M = (sol_K27M[i,:] .+ sol_K27M[i+1,:])./(sol_K27M[i,1] .+ sol_K27M[i+1,1]);

        f_norm2 = plot([
            scatter(x = time_sol, y = norm_K27M, mode = "lines", marker = attr(color = c_K27M[1]), line_width = 6.0, name = "Solution", visible = true),
            scatter(x = time_sol, y = norm_WT, mode = "lines", marker = attr(color = c_WT[1]), line_width = 6.0, name = "Solution", visible = true)
        ])
    else

        norm_WT = sol_WT[i,:]./sol_WT[i,1];
        norm_K27M = sol_K27M[i+1,:]./sol_K27M[i+1,1]

        f_norm2 = plot([
            scatter(x = time_sol, y = norm_K27M, mode = "lines", marker = attr(color = c_K27M[1]), line_width = 6.0, name = "Solution", visible = true),
            scatter(x = time_sol, y = norm_WT, mode = "lines", marker = attr(color = c_WT[1]), line_width = 6.0, name = "Solution", visible = true)
         ])
    end

    ymin = minimum(vcat(norm_WT,norm_K27M));
    ymax = maximum(vcat(norm_WT,norm_K27M)); 
    ymin_pad, ymax_pad, y_ticks, y_step = nice_ticks(ymin, ymax, n_ticks = 5);
        
    width_f=375; height_f=375;

    index = findfirst(==(namePops_WT[i]), namePops_WT);

    if index <= 3
        xaxis_attr = attr(
            showticklabels = false,
            tickvals = [28, 56, 84, 112],
            range = [20, 120],
            tickfont = attr(size = 28)
        );
        yaxis_attr = attr(
            showticklabels = true,
            range = [ymin_pad, ymax_pad],
            tickmode = "array",
            tickvals = y_ticks,
            tickfont = attr(size = 28)
        );
    else
        xaxis_attr = attr(
            showticklabels = true,
            tickvals = [28, 56, 84, 112],
            range = [20, 120],
            tickfont = attr(size = 28)
        );
        yaxis_attr = attr(
            showticklabels = true,
            range = [ymin_pad, ymax_pad],
            tickmode = "array",
            tickvals = y_ticks,
            tickfont = attr(size = 28)
        );
    end

    f_norm2.plot.layout.xaxis = xaxis_attr;
    f_norm2.plot.layout.yaxis = yaxis_attr;

    relayout!(f_norm2, 
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

    display(f_norm2)
    # figpath = pwd()*"/";
    # savefig(f_norm2,figpath*namePops_WT[i]*" - Normalized t0.svg";width=width_f, height=height_f);
end

## =============================
# Figure S5
# ==============================

pbestu_WT = pbest_WT[idx_unique_WT];
pbestu_K27M = pbest_K27M[idx_unique_K27M];

for i in eachindex(uni_WT)
    xmin = maximum(vcat(0.0,minimum(unique_WT[:,i])));
    xmax = maximum(unique_WT[:,i]);

    xmin_pad, xmax_pad, ticks, step = nice_ticks(xmin, xmax);

    # determine decimal formatting for x axis
    decimals = tick_decimals(step);
    tick_text = tick_label.(ticks, decimals);

    # histogram binning
    nbins = 15;
    bin_size = (xmax_pad - xmin_pad) / nbins;

    yvals = histogram_y_extent(unique_WT[:,i], xmin_pad, xmax_pad, bin_size);
    ymax = maximum(yvals);

    # Histogram
    his = plot(
        histogram(
            x = unique_WT[:,i],
            xbins = attr(
                start = xmin_pad,
                stop  = xmax_pad,
                size  = bin_size
            ),
            showlegend = false,
            marker = attr(
                color = c_WT[1],
                opacity = 0.30,
                line = attr(width=2.0, color=c_WT[1])
            )
        )
    );

    y_min_pad, y_max_pad, y_ticks, y_step = nice_ticks(0.0, ymax);
    y_decimals = tick_decimals(y_step);
    y_tick_text = tick_label.(y_ticks, y_decimals);

    xaxis_attr = attr(
        range = [xmin_pad, xmax_pad],
        tickmode = "array",
        tickvals = ticks,
        ticktext = tick_text,
        tickformat = "",
        automargin = true,
        tickfont = attr(size = 28),
        showgrid = false
    );

    yaxis_attr = attr(
        range = [y_min_pad, y_max_pad],
        tickmode = "array",
        tickvals = y_ticks,
        ticktext = y_tick_text,
        tickformat = "",
        tickfont = attr(size = 28)
    );

    his.plot.layout.xaxis = xaxis_attr;
    his.plot.layout.yaxis = yaxis_attr;

    width_f=375; height_f=375;

    annotations =[
        attr(
            text = uni_WT[i],
            x = 0.5,
            y = 1.0,
            xref = "paper",
            yref = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = false,
            font = attr(size = 32)
        )
    ];

    relayout!(
        his,
        template = craig_lab_template,
        annotations = annotations,
        showlegend = false,
        margin = attr(l=55, r=15, t=45, b=45),
        width = width_f,
        height = height_f,
        barmode = "overlay",
        autosize = false,
        shapes = [
            vline(
                pbestu_WT[i],
                line = attr(
                    color = c_WT[2],
                    width = 4.0,
                    dash = "dash"
                )
            )
        ]
    );

    display(his)
    # figpath = pwd()*"/";
    # savefig(his,figpath*"WT Param $i.svg";width=width_f, height=height_f);
end

for i in eachindex(uni_K27M)
    xmin = maximum(vcat(0.0,minimum(unique_K27M[:,i])));
    xmax = maximum(unique_K27M[:,i]);

    xmin_pad, xmax_pad, ticks, step = nice_ticks(xmin, xmax);

    # determine decimal formatting for x axis
    decimals = tick_decimals(step);
    tick_text = tick_label.(ticks, decimals);

    # histogram binning
    nbins = 15;
    bin_size = (xmax_pad - xmin_pad) / nbins;

    yvals = histogram_y_extent(unique_K27M[:,i], xmin_pad, xmax_pad, bin_size);
    ymax = maximum(yvals);

    # Histogram
    his = plot(
        histogram(
            x = unique_K27M[:,i],
            xbins = attr(
                start = xmin_pad,
                stop  = xmax_pad,
                size  = bin_size
            ),
            showlegend = false,
            marker = attr(
                color = c_K27M[1],
                opacity = 0.30,
                line = attr(width=2.0, color=c_K27M[1])
            )
        )
    );

    y_min_pad, y_max_pad, y_ticks, y_step = nice_ticks(0.0, ymax);
    y_decimals = tick_decimals(y_step);
    y_tick_text = tick_label.(y_ticks, y_decimals);

    xaxis_attr = attr(
        range = [xmin_pad, xmax_pad],
        tickmode = "array",
        tickvals = ticks,
        ticktext = tick_text,
        tickformat = "",
        automargin = true,
        tickfont = attr(size = 28),
        showgrid = false
    );

    yaxis_attr = attr(
        range = [y_min_pad, y_max_pad],
        tickmode = "array",
        tickvals = y_ticks,
        ticktext = y_tick_text,
        tickformat = "",
        tickfont = attr(size = 28)
    );

    his.plot.layout.xaxis = xaxis_attr;
    his.plot.layout.yaxis = yaxis_attr;

    width_f=375; height_f=375;

    annotations =[
        attr(
            text = uni_K27M[i],
            x = 0.5,
            y = 1.0,
            xref = "paper",
            yref = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = false,
            font = attr(size = 32)
        )
    ];

    relayout!(
        his,
        template = craig_lab_template,
        annotations = annotations,
        showlegend = false,
        margin = attr(l=55, r=15, t=45, b=45),
        width = width_f,
        height = height_f,
        barmode = "overlay",
        autosize = false,
        shapes = [
            vline(
                pbestu_K27M[i],
                line = attr(
                    color = c_K27M[2],
                    width = 4.0,
                    dash = "dash"
                )
            )
        ]
    );

    display(his)
    # figpath = pwd()*"/";
    # savefig(his,figpath*"K27M Param $i.svg";width=width_f, height=height_f);
end
