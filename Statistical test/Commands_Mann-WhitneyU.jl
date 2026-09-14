using StatsBase, HypothesisTests, CSV, DataFrames

# H3-WT
bloodModel_WT = "Combined HSC & MPP Self-Renewal"; 
include(joinpath(@__DIR__, "..", "Models", bloodModel_WT * ".jl"));
df_WT = CSV.read(joinpath(@__DIR__,"..", "Results", "Accepted Parameters (WT).csv"), DataFrame);
p_WT = Matrix(df_WT[:,3:end]);
labs_WT = labs;

# H3-K27M
bloodModel_K27M = "Transiting HSCs & MPP Self-Renewal"; 
include(joinpath(@__DIR__, "..", "Models", bloodModel_K27M * ".jl"));
df_K27M = CSV.read(joinpath(@__DIR__,"..", "Results", "Accepted Parameters (K27M).csv"), DataFrame);
p_K27M = Matrix(df_K27M[:,3:end]);
labs_K27M = labs; 

# ==============================
# Mann-Whitney U test
# ==============================

common = intersect(labs_K27M, labs_WT);
common_order = [
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
index = Dict(l => i for (i, l) in enumerate(common_order));
perm = sortperm(common, by = x -> get(index, x, Inf));
common = common[perm];

idx_WT   = [findfirst(==(x), labs_WT)   for x in common];
idx_K27M = [findfirst(==(x), labs_K27M) for x in common];

common_WT   = p_WT[:,idx_WT];
common_K27M = p_K27M[:,idx_K27M];

results = DataFrame(
    parameter = String[],
    effect_size = Float64[],     
    rank_biserial = Float64[],
    pvalue = Float64[]
);

for j in eachindex(common)
    x = common_WT[:,j]; y = common_K27M[:,j];

    # Compare K27M to wild-type
    MWU = MannWhitneyUTest(y,x);
    pval = pvalue(MWU);

    ranks = ordinalrank(vcat(y,x));
    n1 = length(y); n2 = length(x);

    R1 = sum(ranks[1:n1]);
    R2 = sum(ranks[n1+1:n1+n2]);

    U1 = R1 - n1*(n1+1)/2;
    U2 = R2 - n2*(n2+1)/2;

    #Common language effect size
    f = U1/(n1*n2);

    # Rank-biserial correlation
    r = 2*f - 1;

    push!(results, (
        parameter = common[j],
        effect_size = f,
        rank_biserial = r,
        pvalue = pval
    ))
end

# CSV.write(joinpath(@__DIR__, "..", "Results", "Mann-Whitney U test.csv"),results);
