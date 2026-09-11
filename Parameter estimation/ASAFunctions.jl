function inbounds(x, L, U)
    for i in eachindex(x)
        if x[i] < L[i] || x[i] > U[i]
            return false
        end
    end
    return true
end

function ASA(theta0,L,U,prob,costdata::CostData; rng=Random.default_rng())
    #Adaptive simulated annealing algorithm created by Ingber (1989) and adapted from Chen & Luk (1999)

    # ---------------------------------------------------------------------- %
    # Initializing 
    
    # Create integrator for ODE solve
    integrator = init(prob, Rodas5(); saveat=costdata.save_times, abstol=1e-8, reltol=1e-6);
    # init(prob, AutoTsit5(Rosenbrock23()); saveat=costdata.save_times, abstol=1e-8, reltol=1e-6);

    # Cost function of initial point
    theta_best = copy(theta0);
    CF_best = costfunction!(integrator, theta_best, costdata);

    theta_new = copy(theta_best);
    CF_new = CF_best;
    theta_pert = copy(theta_best);
    
    # Fixed parameter values
    N  = length(theta0);  # Parameter space dimension
    m  = 5;               # Reannealing parameter

    Ta_0 = CF_best;       # Acceptance temperature at na = 0
    Ta   = Ta_0;          # Acceptance temperature
    T_0  = ones(N);       # Generating temperature at n = 0
    T    = T_0;           # Generating temperature
    na   = 0;             # Acceptance annealing time
    n    = zeros(N);      # Generating annealing time  
    eps  = 0.005;

    # Initializing algorithm variables
    N_gen = 1;
    N_acc = 1;
    N_acc_total = 1;
    e_val = 1;

    reannealing_steps = 0;
    CF_test = Vector{Float64}();

    max_eval = 500000;
    escape   = 0;           # Escape value (to stop ASA)
    FT       = 0.05;        # Function tolerance

# ---------------------------------------------------------------------- %
    # ASA Algorithm 
    
    while escape < 1 && e_val < max_eval

        # Generate a new point
        u = rand(rng, Uniform(0,1),N);
        g = sign.(u.-0.5).*T.*((inv.(T).+1).^(abs.(2*u.-1)).-1);
        
        theta_gen = theta_best .+ g.*(U - L);
        
        # If needed, repeat until the point is in bounds
        while !inbounds(theta_gen, L, U)
            u = rand(rng, Uniform(0,1),N);
            g = sign.(u.-0.5).*T.*((inv.(T).+1).^(abs.(2*u.-1)).-1)
            
            theta_gen = theta_best .+ g.*(U - L);
        end
        
        # Update number of generated points
        N_gen += 1;
        
        # Evaluate cost function of the generated point
        CF_gen = costfunction!(integrator, theta_gen, costdata);
        e_val += 1;
        
        if CF_gen < CF_best
            # accept the new point
            CF_new = CF_gen;
            theta_new .= theta_gen;

            # set point as best value
            CF_best = CF_new;
            theta_best .= theta_new;
            #println(CF_best)
    
            # Update number of accepted points
            N_acc += 1;
            N_acc_total += 1;
    
        else
            q = rand(rng, Uniform(0,1));
            Pa = exp(-(CF_gen - CF_best)/Ta); # standand Metropolis acceptance (less conservative)
            #Pa = 1/(1 + exp((CF_gen - CF_best)/Ta)) # logistic acceptance (more conservative)
    
            if q <= Pa
                # accept the new point
                CF_new = CF_gen;
                theta_new .= theta_gen;
    
                # Update number of accepted points
                N_acc += 1;
                N_acc_total += 1;
            end
        end
    
        if N_acc >= 30    
            # Reannealing takes place
            # println("Reannealing!")

            s = zeros(N);
            theta_pert .= theta_best;
            
            for k = 1:N
                val = theta_pert[k];
                delta = eps * max(abs(val), 1e-8);
                theta_pert[k] = val + delta # relative perturbation (sensitivity is more consistently measured across parameters)
                # theta_pert[k] = val + eps; # finite perturbation (bad for parameters that differ in magnitude)
                
                CF_pert = costfunction!(integrator, theta_pert, costdata);
                s[k] = abs((CF_pert - CF_best)/delta);
                # s[k] = abs((CF_pert - CF_best)/eps);

                # reset best parameter value
                theta_pert[k] = val;
            end
            s_max = maximum(s);
            
            T_0 = ones(N);
            for k = 1:N
                T[k] = (s_max/s[k]).*T[k];
            end
            n = (-(1/m).*log.(T./T_0)).^N;
            zero_idx = findall(idx -> idx < 0 ,n);
            n[zero_idx] .= 0

            Ta_0 = CF_new;   # set as last accepted point
            Ta   = CF_best;  # set as best point
            na   = (-(1/m)*log(Ta/Ta_0)).^N;
            
            # Reinitializing after reannealing
            N_acc = 1;
    
            # Checking if the algorithm stops after reannealing
            reannealing_steps += 1;
            push!(CF_test, CF_best);

            if reannealing_steps >= 5 
                # println("Change in CF in last 5 reannealing steps: ", CF_test[reannealing_steps-4] - CF_test[reannealing_steps])
                if CF_test[reannealing_steps-4] - CF_test[reannealing_steps] <= FT
                    escape = 1;
                else
                    escape = 0;
                end
            end
    
        end 
        
        if N_gen >= 800
            # Temperature annealing takes place

            n  = n .+ 1;
            T  = T_0.*exp.(-m.*n.^(1/N));
            na = na + 1;
            Ta = Ta_0*exp(-m*na^(1/N));
    
            # Reinitializing after annealing
            N_gen = 1;
        end
    end

    # println("Best CF from run: ", CF_best)
    # println("Total cost evaluations: ", e_val)
    # println("Total accepted points: ", N_acc_total)
    # println("Acceptance ratio: ", (N_acc_total/e_val)*100)

    return (theta_best, CF_best)

end

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
