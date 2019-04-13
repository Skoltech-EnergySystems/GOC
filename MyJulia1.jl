# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/src/");
# push!(LOAD_PATH, "/Users/ivan/Documents/PhD/GOC/");
using Revise
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
opf = create_model(PN)
# Writing into the file "solution1.txt"
Bus_I = OPF[1]
V_I = OPF[2]
Theta_I = OPF[3]
BCS_I = OPF[4]
Gen_I = OPF[5]
Gen_ID = OPF[6]
P_G = OPF[7]
Q_G = OPF[8]

NBus = size(Bus_I,1)
sSh_ind = [s.i for s in PN.sShuntList]
NGen = size(Gen_I,1)

Nod = zeros(NBus,4)
Nod[1:NBus,1] = Bus_I[:,1]
Nod[1:NBus,2] = V_I[:,1]
Nod[1:NBus,3] = Theta_I[:,1]
Nod[sSh_ind,4] = BCS_I

Gen = zeros(NGen,4)
Gen[1:NGen,1] =Gen_I[:,1]
Gen[1:NGen,2] =Gen_ID[:,1]
Gen[1:NGen,3] =P_G[:,1]
Gen[1:NGen,4] =Q_G[:,1]

    open("solution1.txt", "w") do f1
        write(f1, "--bus section\r\n")
        write(f1, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\r\n")
        for i in 1:NBus
            write(f1, join(Nod[i,:], ","))
            write(f1, "\r\n")
        end
        write(f1, "--generator section\r\n")
        write(f1, "i, id, p(MW), q(MVAR)\r\n")
        for i in 1:NGen
            write(f1, join(Gen[i,:], ","))
            write(f1, "\r\n")
        end
    end
end
