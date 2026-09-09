# ======================================= 
# Helper functions for plotting figures
# =======================================

# ---------------------------------------------------------------------- %
function hex2rgba(hex::String, opacity::Float64)
    
    # Remove the "#"
    hex = replace(hex, "#" => "")
    
    # Parse the Red, Green, and Blue pairs
    r = parse(Int, hex[1:2], base=16)
    g = parse(Int, hex[3:4], base=16)
    b = parse(Int, hex[5:6], base=16)
    
    # Return rgba stringH3
    return "rgba($r, $g, $b, $opacity)"
end

# ---------------------------------------------------------------------- %
function GetBestParam(p,CF)

    CF_best = minimum(CF);
    idx = getindex.(indexin(CF_best,CF),1);
    p_best = p[idx,:];

    bestfit = [p_best, CF_best];
    return bestfit
end

# ---------------------------------------------------------------------- %
function SolveBloodODE(p_best,time,IC)

    tspan = (minimum(time), maximum(time));
    stime = minimum(time):0.5:maximum(time);

    prob = ODEProblem(bloodODE!,IC,tspan,p_best);
    sol = solve(prob, Tsit5(), dt = 0.5, saveat = stime);
    
    return sol
end

# ---------------------------------------------------------------------- %
function nice_ticks(xmin, xmax; n_ticks::Integer = 4)

    xmax = max(xmax, xmin + 1e-12)
    span = xmax - xmin
    steps = n_ticks - 1

    nice_steps = Float64[]

    for exp in -3:6
        append!(nice_steps, [1.0, 1.25, 2.0, 2.5, 4.0, 5.0, 6.0] .* 10.0^exp)
    end

    sort!(nice_steps)

    # Smallest step giving exactly n_ticks and covers data with padding
    step = nothing
    xmin_pad = 0.0

    for candidate in nice_steps
        if steps * candidate + 1e-12 < span
            continue
        end

        # candidate range can cover the data if we align it with the minimum
        left_start = max(0.0, floor(xmin / candidate) * candidate)
        if left_start + steps * candidate >= xmax - 1e-12
            step = candidate
            xmin_pad = left_start
            break
        end

        # or shift the range forward to just cover the maximum while still starting before the minimum
        right_start = max(0.0, ceil((xmax - steps * candidate) / candidate) * candidate)
        if right_start <= xmin + 1e-12
            step = candidate
            xmin_pad = right_start
            break
        end
    end

    if isnothing(step)
        step = nice_steps[end]
        xmin_pad = max(0.0, floor(xmin / step) * step)
    end

    xmax_pad = xmin_pad + steps * step

    # exactly n_ticks evenly spaced tick marks
    ticks = xmin_pad .+ step .* (0:steps)

    return xmin_pad, xmax_pad, ticks, step
end

# ---------------------------------------------------------------------- %
function tick_decimals(step)
    if step >= 1.0 || step <= 0.0
        return 0
    end
    decimals = 0
    x = step
    while !isapprox(round(x), x; atol=1e-14, rtol=1e-14) && decimals < 12
        x *= 10
        decimals += 1
    end
    return decimals
end

# ---------------------------------------------------------------------- %
function tick_label(x, decimals)
    if decimals == 0
        return string(Int(round(x)))
    end
    s = string(round(x, digits=decimals))
    if endswith(s, ".0")
        return s[1:end-2]
    end
    return s
end

# ---------------------------------------------------------------------- %
function histogram_y_extent(x, start, stop, size; histnorm="")
    n_bins = Int(round((stop - start) / size))
    n_bins = max(n_bins, 1)
    counts = zeros(Float64, n_bins)

    for value in x
        if value < start || value > stop
            continue
        end
        idx = Int(floor((value - start) / size)) + 1
        idx = clamp(idx, 1, n_bins)
        counts[idx] += 1
    end

    total = length(x)
    if histnorm == "percent"
        return counts ./ total .* 100
    elseif histnorm == "density" || histnorm == "probability density"
        return counts ./ (total * size)
    elseif histnorm == "probability"
        return counts ./ total
    else
        return counts
    end
end

# ---------------------------------------------------------------------- %
function summary_trace(x, mean, se; width=0.45)
    [
        # vertical error bar
        scatter(
            x = [x, x],
            y = [mean-se, mean+se],
            mode = "lines",
            line = attr(color="black", dash = "dot", width=2),
            hoverinfo = "skip",
            showlegend = false
        ),

        # horizontal mean line
        scatter(
            x = [x-width/2, x+width/2],
            y = [mean, mean],
            mode = "lines",
            line = attr(color="black", width=3),
            hoverinfo = "skip",
            showlegend = false
        ),

        # error-bar caps
        scatter(
            x = [x-width/4, x+width/4],
            y = [mean+se, mean+se],
            mode = "lines",
            line = attr(color="black", width=2),
            hoverinfo = "skip",
            showlegend = false
        ),

        scatter(
            x = [x-width/4, x+width/4],
            y = [mean-se, mean-se],
            mode = "lines",
            line = attr(color="black", width=2),
            hoverinfo = "skip",
            showlegend = false
        )
    ]
end
