include("parsing.jl")
include("structure.jl")

rawFile = "./data/case.raw"
contFile = "./data/case.con"
costsFile = "./data/case.rop"
inlFile = "./data/case.inl"


println("Starting pre-processing.\n Parsing Data and storing it in DataFrames")
println("------------------------------")
# process raw Data
println("raw data at $rawFile")
Data = initialize_data_struct()
@time raw_parser!(rawFile, Data);
println("------------------------------")
# process contingencies
println("file with contingencies at $contFile")
Contin = initialize_contin_struct()
@time cont_parser!(contFile, Contin)
println("------------------------------")
# process Costs
println("file with costs at $contFile")
Costs = initialize_costs_struct()
@time costs_parser!(costsFile, Costs)
println("------------------------------")
# process *.inl data
println("case.inl at $inlFile")
Resp = initialize_response_struct()
response_parser!(inlFile, Resp)
println("------------------------------")
