# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/src/");
# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/");


cur_path = pwd();
code_path = cur_path * "/src/";
data_path = cur_path * "/data/";

using Revise
includet(code_path * "Parsing.jl")
includet(code_path * "DataStructures.jl")
includet(code_path * "NetworkData.jl")

Revise.track(code_path * "Parsing.jl")
Revise.track(code_path * "DataStructures.jl")
Revise.track(code_path * "NetworkData.jl")

# rawFile = data_path * "case.raw"
# contFile = data_path * "case.con"
# costsFile = data_path * "case.rop"
# inlFile = data_path * "case.inl"
data_pathes = Dict(:raw => data_path * "case.raw",
                    :contin => data_path * "case.con",
                    :costs => data_path * "case.rop",
                    :inl => data_path * "case.inl"
)

PN, mainData, costsData, continData = parser(data_pathes)


# include("model.jl")
# opf = create_model(PN)
