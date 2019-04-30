# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/src/");
# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/");

using Revise

# paths manipulations
# Julia 1.1.0
#=
Path = splitpath(pwd())
i = findlast(s -> s == "GOC", Path)
proj_path = joinpath(Path[1:i]...)
=#
# Julia 1.0 and older
cur_path, goc = splitdir(pwd())
while goc != "GOC"
    cur_path, goc = splitdir(cur_path)
end
proj_path = joinpath(cur_path, goc);
#
code_path = joinpath(proj_path, "src");
data_path = joinpath(proj_path, "data");

# Revise tracking
includet(joinpath(code_path , "Parsing.jl"));
includet(joinpath(code_path , "DataStructures.jl"));
includet(joinpath(code_path , "NetworkData.jl"));
includet(joinpath(code_path , "model.jl"));
Revise.track(joinpath(code_path , "Parsing.jl"));
Revise.track(joinpath(code_path , "DataStructures.jl"));
Revise.track(joinpath(code_path , "NetworkData.jl"));
Revise.track(joinpath(code_path , "model.jl"));


data_paths = Dict(:raw => joinpath(data_path , "case.raw"),
                    :contin => joinpath(data_path , "case.con"),
                    :costs => joinpath(data_path , "case.rop"),
                    :inl => joinpath(data_path , "case.inl")
);

# parsing
PN, mainData, costsData, continData = parser(data_paths);


# include("model.jl")
# opf = create_model(PN)

shunt_idx = [S.i for S in PN.sShuntList ]
