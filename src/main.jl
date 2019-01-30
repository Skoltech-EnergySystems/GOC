include("Parsing.jl")
include("DataStructures.jl")
include("NetworkData.jl")
include("network_init.jl")

rawFile = "./data/case.raw"
contFile = "./data/case.con"
costsFile = "./data/case.rop"
inlFile = "./data/case.inl"


println("Starting pre-processing.\n Parsing Data and storing it in DataFrames")
println("------------------------------")
# process raw Data
println("raw data at $rawFile")
mainData = initialize_data_struct()
@time raw_parser!(rawFile, mainData);
println("------------------------------")
# process contingencies
println("file with contingencies at $contFile")
continData = initialize_contin_struct()
@time cont_parser!(contFile, continData)
println("------------------------------")
# process Costs
println("file with costs at $contFile")
costsData = initialize_costs_struct()
@time costs_parser!(costsFile, costsData)
println("------------------------------")
# process *.inl data
# println("case.inl at $inlFile")
# respData = initialize_response_struct()
# @time response_parser!(inlFile, respData)
# println("------------------------------")

# PNetwork initialization with Data
PN = PNetwork()
@time network_init(PN, mainData)

# PCosts init
include("network_init.jl")
PC = PCosts()
costs_init(costsData, PC, PN.s)
