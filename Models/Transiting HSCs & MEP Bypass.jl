# Initializing model fit 

include(@__DIR__, "..", "Data","BloodData.jl");

cellCount_WT   = [HSC1_WT HSC2_WT MPP_WT CMP_WT GMP_WT MEP_WT EPp_WT EPn_WT MLP_WT];
cellCount_K27M = [HSC1_K27M HSC2_K27M MPP_K27M CMP_K27M GMP_K27M MEP_K27M EPp_K27M EPn_K27M MLP_K27M];

namePops = ["HSC1","HSC2","MPP","CMP","GMP","MEP","EP CD71+","EP CD71-","MLP"];
nPops = size(cellCount_WT,2);

N_WT   = size(cellCount_WT,1)*size(cellCount_WT,2);
N_K27M = size(cellCount_K27M,1)*size(cellCount_K27M,2);

labs = ["a<sub>HSC1</sub>", "p<sub>HSC1\u2192MPP</sub>", "a<sub>HSC2</sub>", "p<sub>HSC2\u2192MPP</sub>", "p<sub>MPP\u2192CMP</sub>", "p<sub>CMP\u2192GMP</sub>", "p<sub>CMP\u2192MEP</sub>", "p<sub>MPP\u2192MLP</sub>", 
"p<sub>MEP\u2192EP+</sub>", "p<sub>EP+\u2192EP-</sub>", "p<sub>MPP\u2192EP+</sub>", "d<sub>MPP</sub>", "d<sub>CMP</sub>", "p<sub>GMP\u2192M</sub>", "d<sub>MEP</sub>", "p<sub>MLP\u2192M</sub>", "d<sub>EP+</sub>", "p<sub>EP-\u2192M</sub>", 
"r<sub>HSC1\u2192HSC2</sub>", "r<sub>HSC2\u2192HSC1</sub>"];

N_param = size(labs,1);

dCell_WT = Vector{Matrix{Any}}(undef, nPops);
dCell_K27M = Vector{Matrix{Any}}(undef, nPops);

i28 = findall(==(28), time_d); i56 = findall(==(56), time_d); 
i84 = findall(==(84), time_d); i112 = findall(==(112), time_d); 
for i = 1:nPops
    dCell_WT[i] = permutedims(
        hcat(cellCount_WT[:,i][i28], cellCount_WT[:,i][i56], cellCount_WT[:,i][i84], cellCount_WT[:,i][i112])
        );
    dCell_K27M[i] = permutedims(
        hcat(cellCount_K27M[:,i][i28], cellCount_K27M[:,i][i56], cellCount_K27M[:,i][i84], cellCount_K27M[:,i][i112])
    );
end

# ODE initial condition
avgCell_WT = zeros(size(mean_WT,1),nPops);
avgCell_K27M = zeros(size(mean_K27M,1),nPops);
stdCell_WT = zeros(size(mean_WT,1),nPops);
stdCell_K27M = zeros(size(mean_K27M,1),nPops);
IC_WT = zeros(1,nPops);
IC_K27M = zeros(1,nPops);

for i = 1:nPops 
    avgCell_WT[:,i]   = mean(dCell_WT[i],dims=2);
    avgCell_K27M[:,i] = mean(dCell_K27M[i],dims=2);
    stdCell_WT[:,i]   = std(dCell_WT[i],dims=2);
    stdCell_K27M[:,i] = std(dCell_K27M[i],dims=2);

    IC_WT[1,i] = avgCell_WT[1,i];
    IC_K27M[1,i] = avgCell_K27M[1,i];
end


# Generating initial points in parameter space
#            [aHSC1  p1_MPP aHSC2  p2_MPP pM_CMP pC_GMP pC_MEP pM_MLP pME_EP pEP_EP pM_EPp dMPP   dCMP   dGMP   dMEP   dMLP   dEP+   dEP-   r12    r21   ];
L_WT   = vec([0.0000 1.0000 0.0000 1.5000 1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000]);
U_WT   = vec([0.0200 8.0000 0.0400 8.0000 8.0000 4.0000 5.0000 0.2000 10.000 3.0000 5.0000 4.0000 2.0000 1.0000 4.0000 0.1000 2.0000 2.0000 0.0500 0.0200]);

L_K27M = vec([0.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000]);
U_K27M = vec([0.0500 7.0000 0.0500 8.0000 8.0000 4.0000 5.0000 0.1000 10.000 1.5000 6.0000 4.0000 2.0000 4.0000 4.0000 0.2000 2.0000 4.0000 0.0500 0.0500]);
#-----------------------------------------------------------------------------------

function bloodODE!(dy,y,param,t)
# Model where we consider HSC1 and HSC2 as transiting populations of LT-HSCs.
    
        a1  = param[1];  p1  = param[2];
        a2  = param[3];  p2  = param[4]
        p3  = param[5];  p4  = param[6];
        p5  = param[7];  p6  = param[8];
        p7  = param[9];  p8  = param[10];
        p9  = param[11]; d3  = param[12];
        d4  = param[13]; d5  = param[14];
        d6  = param[15]; d7  = param[16];
        d8  = param[17]; d9  = param[18];
        r12 = param[19]; r21 = param[20];
        
        HSC1 = y[1];
        HSC2 = y[2]; MPP = y[3];     
        CMP  = y[4]; GMP = y[5];     
        MEP  = y[6]; EPp = y[7]; 
        EPn  = y[8]; MLP = y[9];
        
        dHSC1 = a1*HSC1 + r21*HSC2;
        dHSC2 = a2*HSC2 + r12*HSC1;
        dMPP  = p1*HSC1 + p2*HSC2 - (d3+p3+p6+p9)*MPP;
        dCMP  = p3*MPP - (d4+p4+p5)*CMP;
        dGMP  = p4*CMP - d5*GMP;
        dMEP  = p5*CMP - (d6+p7)*MEP;
        dEPp  = p7*MEP + p9*MPP - (d8+p8)*EPp;
        dEPn  = p8*EPp - d9*EPn;
        dMLP  = p6*MPP - d7*MLP;
    
        dy[1] = dHSC1
        dy[2] = dHSC2; dy[3] = dMPP;     
        dy[4] = dCMP;  dy[5] = dGMP;     
        dy[6] = dMEP;  dy[7] = dEPp; 
        dy[8] = dEPn;  dy[9] = dMLP;
    
end

struct CostData
    cellCount::Matrix{Float64}
    save_times::Vector{Float64}
    idxs::NTuple{4, Vector{Int}}   # indexes (i28, i56, i84, i112)
end

function costfunction!(integrator, param::AbstractVector, data::CostData)

    cellCount   = data.cellCount;
    save_times  = data.save_times;
    (i28, i56, i84, i112) = data.idxs;

    # update parameters in-place
    integrator.p .= param

    # reset state to t0 and initial condition
    reinit!(integrator)

    solve!(integrator)

    sol = integrator.sol;

    nCellPop = size(cellCount,2);
    NSE = 0.0;

    for i = 1:nCellPop
        data_block  = @view cellCount[:,i]
        model_block = sol[i,:]

        d1 = data_block[i28]
        d2 = data_block[i56]
        d3 = data_block[i84]
        d4 = data_block[i112]

        data_mat = hcat(d1,d2,d3,d4)'
        μ = mean(data_mat)

        NSE += sum(((data_mat .- model_block) ./ μ).^2)
    end

    return isfinite(NSE) ? NSE : Inf
end