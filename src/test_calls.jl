workspace()


include("Parsing.jl")
include("DataStructures.jl")
include("NetworkData.jl")
# include("network_init.jl")

rawFile = "./data/case.raw"
contFile = "./data/case.con"
costsFile = "./data/case.rop"
inlFile = "./data/case.inl"

PN = PNetwork()
mainData = PNetworkData()
@time raw_parser!(rawFile, mainData, PN);

costsData = CostsStruct()
@time costs_parser!(costsFile, costsData, PN)

include("Parsing.jl")
# respData = initialize_response_struct()
response_parser!(inlFile, PN)#, respData)


continData = ContingenciesStruct()
@time cont_parser!(contFile, continData)


include("model.jl")
opf = create_model(PN)
