workspace()

push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/Parsing.jl")


using Revise
Revise.track("Parsing.jl")

includet("Parsing.jl")
includet("DataStructures.jl")
includet("NetworkData.jl")
# include("network_init.jl")

rawFile = "./data/case.raw"
contFile = "./data/case.con"
costsFile = "./data/case.rop"
inlFile = "./data/case.inl"

# create an empty Network structure and fill it with parsed data
PN = PNetwork()
mainData = PNetworkData() # empty structure with DataFrames to store data
@time raw_parser!(rawFile, mainData, PN);

costsData = CostsStruct()
@time costs_parser!(costsFile, costsData, PN)

include("Parsing.jl")
# respData = initialize_response_struct()
response_parser!(inlFile, PN)#, respData)


continData = ContingenciesStruct()
@time cont_parser!(contFile, continData)


# include("model.jl")
# opf = create_model(PN)
