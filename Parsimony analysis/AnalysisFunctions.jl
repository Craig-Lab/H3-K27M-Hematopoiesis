# ---------------------------------------------------------------------- %
# Functions for getting the best fit solution

function GetBestParam(p,CF)

    CF_best = minimum(CF);
    idx = getindex.(indexin(CF_best,CF),1);
    p_best = p[idx,:];

    bestfit = [p_best, CF_best];
    return bestfit
end

function GetBestIdx(p,CF)

    CF_best = minimum(CF);
    idx = getindex.(indexin(CF_best,CF));
    p_best = p[idx,:];

    bestfit = [p_best, CF_best];
    return bestfit
end


function AICc(RSS,k,N)
# Returns the AICc using the residual sum of squares (RSS), the number of parameters (k), and the number of observations (N).
        
    if k < N/40
        AICc = N*log(RSS/N) + 2*k;
    end 
    if k >= N/40
        AICc = N*log(RSS/N) + (2*k*N)/(N-k-1); 
    end
    return AICc
    
end

function BIC(RSS,k,N)
    # Returns the BIC using the residual sum of squares (RSS), the number of parameters (k), and the number of observations (N).
        
    BIC = N*log(RSS/N) + k*log(N);
    return BIC
        
end

function m2LL(RSS,cellCount)
# Returns the log-likelihood using the residual sum of squares (RSS) and the number of observations (N).
    
    Npop = length(cellCount); 
    MVS = mean(cellCount,dims = 1);
    m2LogL = RSS + Npop*sum(log.((2*pi).*MVS))
    return m2LogL
        
end