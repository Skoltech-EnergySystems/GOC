# Function MyJulia1

function MyJulia1(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
        contFile = InFile1
        inlFile = InFile2
        rawFile = InFile3
        costsFile = InFile4

### ====> LOADING FUNCTIONS
println("\n ... LOADING FUNCTIONS")
code_path = joinpath(pwd(), "src")
include(joinpath(code_path , "Parsing.jl"));
include(joinpath(code_path , "DataStructures.jl"));
include(joinpath(code_path , "NetworkData.jl"));
include(joinpath(code_path , "model.jl"))

### ====> PARSING DATA
println("\n ... PARSING DATA")
data_paths = Dict(:raw => rawFile,
                    :contin => contFile,
                    :costs => costsFile,
                    :inl => inlFile);
PN, mainData, costsData, continData = @time parser(data_paths);

### ====> SOLVING OPF
println("\n ... SOLVING OPF")
println("\n ... creating the model")
Bus_I, V_I, Theta_I, BCS_I, Gen_I, Gen_ID, P_G, Q_G = @time create_model(PN)

# println("\n ... solving the model")
# res = @time optimize!(opf)

### ====> Writing solution1.txt
println("\n ... WRITING solution1.txt")
@time  include("write_solution1.jl")

end
