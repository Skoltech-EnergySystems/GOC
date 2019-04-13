# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/src/");
# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/");
# using Revise
function MyJulia1(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
        rawFile = InFile1
        contFile = InFile2
        costsFile = InFile3
        inlFile = InFile4
# paths manipulations
# Julia 1.1.0
#=
Path = splitpath(pwd())
i = findlast(s -> s == "GOC", Path)
proj_path = joinpath(Path[1:i]...)
=#
# Julia 1.0 and older


### ====> LOADING FUNCTIONS
println("\n ... LOADING FUNCTIONS")

cur_path, goc = splitdir(pwd())
while goc != "GOC"
    cur_path, goc = splitdir(cur_path)
end
proj_path = joinpath(cur_path, goc);
#
code_path = joinpath(proj_path, "src");
data_path = joinpath(proj_path, "data");

# Loading
include(joinpath(code_path , "Parsing.jl"));
include(joinpath(code_path , "DataStructures.jl"));
include(joinpath(code_path , "NetworkData.jl"));
include(joinpath(code_path , "model.jl"))
# Revise.track(joinpath(code_path , "Parsing.jl"));
# Revise.track(joinpath(code_path , "DataStructures.jl"));
# Revise.track(joinpath(code_path , "NetworkData.jl"));
# Revise.track(joinpath(code_path , "model.jl"));

##  Available logic CPU
try
    println("CPU processors available:  ", Sys.CPU_THREADS)
    println("Workers active:  ", nworkers())
catch end

data_paths = Dict(:raw => joinpath(data_path , "case.raw"),
                    :contin => joinpath(data_path , "case.con"),
                    :costs => joinpath(data_path , "case.rop"),
                    :inl => joinpath(data_path , "case.inl"));

### ====> PARSING DATA
println("\n ... PARSING DATA")
PN, mainData, costsData, continData = @time parser(data_paths);

### ====> SOLVING OPF
println("\n ... SOLVING OPF")
println("\n ... creating the model")
Bus_I, V_I, Theta_I, BCS_I, Gen_I, Gen_ID, P_G, Q_G = @time create_model(PN)

# println("\n ... solving the model")
# res = @time optimize!(opf)

end
