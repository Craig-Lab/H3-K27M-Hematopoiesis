using PlotlyJS, LaTeXStrings, CSV, DataFrames, DifferentialEquations, OrdinaryDiffEq, Statistics

include("AnalysisFunctions.jl");

bloodModels = ["Combined HSC", "Combined HSC & MPP Self-Renewal", "Combined HSC & MEP Bypass", "Combined HSC & MPP Self-Renewal & MEP Bypass",
            "Transiting HSCs", "Transiting HSCs & MPP Self-Renewal", "Transiting HSCs & MEP Bypass", "Transiting HSCs & MPP Self-Renewal & MEP Bypass"];
   
conds = ["WT","K27M"];

for j in eachindex(conds)

        Stat = zeros(size(bloodModels,1),5);  
        cond = conds[j];

        for i in eachindex(bloodModels)

            include(joinpath(@__DIR__, "..", "Models", bloodModels[i]*".jl"));

            if cond == "WT"
                df_cond = CSV.read(joinpath(@__DIR__, "..", "Results", bloodModels[i], "WT ASA.csv"), DataFrame);
                cC = cellCount_WT; Nobs = N_WT;
            elseif cond == "K27M"
                df_cond = CSV.read(joinpath(@__DIR__, "..", "Results", bloodModels[i], "K27M ASA.csv"), DataFrame);
                cC = cellCount_K27M; Nobs = N_K27M;
            end

            k_cond = N_param; 
            p_cond = Matrix(df_cond[!,4:end]);
            CF_cond = df_cond[!,"CF"];

            #Statistics for all models
            RSS_cond = GetBestParam(p_cond,CF_cond)[2];

            Stat[i,1] = k_cond;
            Stat[i,2] = AICc(RSS_cond,k_cond,Nobs);
            Stat[i,3] = BIC(RSS_cond,k_cond,Nobs);
            Stat[i,4] = m2LL(RSS_cond,cC);        
            Stat[i,5] = RSS_cond;  
        end

        minAICc = minimum(Stat[:,2]);
        minBIC = minimum(Stat[:,3]);

        dAICc = Stat[:,2] .- minAICc;
        dBIC = Stat[:,3] .- minBIC;

        statTable = DataFrame("Model" => bloodModels, 
                            "N" => Int.(Stat[:,1]),"AICc" => dAICc,"BIC" => dBIC,"CF" => Stat[:,5],"-2LL" => Stat[:,4]
                            );

        show(statTable, allrows=true);

        # CSV.write(joinpath(@__DIR__, "..", "Results", "Model Statistics - "*cond*".csv"), statTable);
end