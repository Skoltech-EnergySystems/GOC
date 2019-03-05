# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/src/");
# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/");

using Revise

"""
the main function which starts all the procedures:
   1. pasing + initialization
   2. model creation ...
"""
 function main()

    # paths manipulations
    Path = splitpath(pwd())
    i = findlast(s -> s == "GOC", Path)
    proj_path = joinpath(Path[1:i]...)
    code_path = joinpath(proj_path, "src");
    data_path = joinpath(proj_path, "data");

    # Revise tracking
    includet(joinpath(code_path , "Parsing.jl"));
    includet(joinpath(code_path , "DataStructures.jl"));
    includet(joinpath(code_path , "NetworkData.jl"));
    Revise.track(joinpath(code_path , "Parsing.jl"));
    Revise.track(joinpath(code_path , "DataStructures.jl"));
    Revise.track(joinpath(code_path , "NetworkData.jl"));

    #  paths to data files
    # TO DO
    # modify to loop through scenario directories
    data_paths = Dict(:raw => joinpath(data_path , "case.raw"),
                        :contin => joinpath(data_path , "case.con"),
                        :costs => joinpath(data_path , "case.rop"),
                        :inl => joinpath(data_path , "case.inl")
    );

    # parsing
    PN, mainData, costsData, continData = parser(data_paths);


    # include("model.jl")
    # opf = create_model(PN)

    return PN, mainData, costsData, continData
end
