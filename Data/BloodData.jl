# Initializing the data
data_WT = dropmissing!(CSV.read("Data - WT.csv", DataFrame));
data_K27M = dropmissing!(CSV.read("Data - K27M.csv", DataFrame));

nExp_WT = count(name -> startswith(String(name), "Exp"), names(data_WT))
nExp_K27M = count(name -> startswith(String(name), "Exp"), names(data_K27M)) 

HSC1_WT = vec(Matrix(data_WT[data_WT.Population .== "HSC1",3:end]));
HSC2_WT = vec(Matrix(data_WT[data_WT.Population .== "HSC2",3:end]));
MPP_WT  = vec(Matrix(data_WT[data_WT.Population .== "MPP",3:end]));
CMP_WT  = vec(Matrix(data_WT[data_WT.Population .== "CMP",3:end]));
GMP_WT  = vec(Matrix(data_WT[data_WT.Population .== "GMP",3:end]));
MEP_WT  = vec(Matrix(data_WT[data_WT.Population .== "MEP",3:end]));
EPp_WT  = vec(Matrix(data_WT[data_WT.Population .== "EP CD71+",3:end]));
EPn_WT  = vec(Matrix(data_WT[data_WT.Population .== "EP CD71-",3:end]));
MLP_WT  = vec(Matrix(data_WT[data_WT.Population .== "MLP",3:end]));

HSC1_K27M = vec(Matrix(data_K27M[data_K27M.Population .== "HSC1",3:end]));
HSC2_K27M = vec(Matrix(data_K27M[data_K27M.Population .== "HSC2",3:end]));
MPP_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "MPP",3:end]));
CMP_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "CMP",3:end]));
GMP_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "GMP",3:end]));
MEP_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "MEP",3:end]));
EPp_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "EP CD71+",3:end]));
EPn_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "EP CD71-",3:end]));
MLP_K27M  = vec(Matrix(data_K27M[data_K27M.Population .== "MLP",3:end]));

time_d = vec(repeat(data_K27M[data_K27M.Population .== "HSC1",2],nExp_K27M));

mean_WT = [
    transpose(mean(Array(data_WT[data_WT.Day .== 28,3:end]), dims=2));
    transpose(mean(Array(data_WT[data_WT.Day .== 56,3:end]), dims=2));
    transpose(mean(Array(data_WT[data_WT.Day .== 84,3:end]), dims=2));
    transpose(mean(Array(data_WT[data_WT.Day .== 112,3:end]), dims=2));;
];
std_WT = [
    transpose(std(Array(data_WT[data_WT.Day .== 28,3:end]), dims=2));
    transpose(std(Array(data_WT[data_WT.Day .== 56,3:end]), dims=2));
    transpose(std(Array(data_WT[data_WT.Day .== 84,3:end]), dims=2));
    transpose(std(Array(data_WT[data_WT.Day .== 112,3:end]), dims=2));;
];
mean_K27M = [
    transpose(mean(Array(data_K27M[data_K27M.Day .== 28,3:end]), dims=2));
    transpose(mean(Array(data_K27M[data_K27M.Day .== 56,3:end]), dims=2));
    transpose(mean(Array(data_K27M[data_K27M.Day .== 84,3:end]), dims=2));
    transpose(mean(Array(data_K27M[data_K27M.Day .== 112,3:end]), dims=2));;
];
std_K27M = [
    transpose(std(Array(data_K27M[data_K27M.Day .== 28,3:end]), dims=2));
    transpose(std(Array(data_K27M[data_K27M.Day .== 56,3:end]), dims=2));
    transpose(std(Array(data_K27M[data_K27M.Day .== 84,3:end]), dims=2));
    transpose(std(Array(data_K27M[data_K27M.Day .== 112,3:end]), dims=2));;
];
mean_time = data_K27M[data_K27M.Population .== "HSC1",2];

