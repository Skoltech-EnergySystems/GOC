### Script -
cd("/Users/david/Dropbox/GitHub/GOC")

###### GENERAL parameters
TimeLimitInSeconds = 600
ScoringMethod = 1
NetworkModel = "Network_01R-10"
sc = 1   # scenario

#########  SIMULATING THE FILE SELECTION
cur_path, goc = splitdir(pwd())
while goc != "GOC"
    cur_path, goc = splitdir(cur_path)
end
proj_path = joinpath(cur_path, goc);
# code_path = joinpath(proj_path, "src");
data_path = joinpath(proj_path, "data", NetworkModel);

contFile =  joinpath(data_path , "scenario_$sc", "case.con")
inlFile =  joinpath(data_path , "case.inl")
rawFile =  joinpath(data_path , "scenario_$sc", "case.raw")
costsFile = joinpath(data_path , "case.rop")


#########  RUNNING as in THE COMPETITION

## First part
include("MyJulia1.jl")
MyJulia1(contFile, inlFile, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)s
