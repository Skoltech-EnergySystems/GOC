# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/src/");
# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/");

using Revise

 function main()

    Path = splitpath(pwd())
    i = findlast(s -> s == "GOC", Path)
    proj_path = joinpath(Path[1:i]...)

    # cur_path = pwd();
    # println(cur_path)
    code_path = joinpath(proj_path, "src");
    data_path = joinpath(proj_path, "data");

    # using Revise
    includet(joinpath(code_path , "Parsing.jl"));
    includet(joinpath(code_path , "DataStructures.jl"));
    includet(joinpath(code_path , "NetworkData.jl"));

    Revise.track(joinpath(code_path , "Parsing.jl"));
    Revise.track(joinpath(code_path , "DataStructures.jl"));
    Revise.track(joinpath(code_path , "NetworkData.jl"));

    # rawFile = data_path * "case.raw"
    # contFile = data_path * "case.con"
    # costsFile = data_path * "case.rop"
    # inlFile = data_path * "case.inl"
    data_paths = Dict(:raw => joinpath(data_path , "case.raw"),
                        :contin => joinpath(data_path , "case.con"),
                        :costs => joinpath(data_path , "case.rop"),
                        :inl => joinpath(data_path , "case.inl")
    );

    PN, mainData, costsData, continData = parser(data_paths);


    # include("model.jl")
    # opf = create_model(PN)
end
