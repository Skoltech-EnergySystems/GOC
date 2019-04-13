### Script
# include("MyJulia1.jl");
#cd("C:/Users/Ильгиз/Documents/Документы_Ильгиз/Skoltech/Grid Competition/Code/julia-knitro-Phase-0-Arun")

# initial model (small)
#=
rawFile = "IEEE14-1_powersystem.raw";
contFile = "IEEE14-1_contingency.csv";
genFile = "IEEE14-1_generator.csv";

MyJulia1(rawFile, genFile, contFile)
=#
#=
# 14 nodes but 72 contingencies
  rawFile = "IEEE14-72_powersystem.raw"
  contFile = "IEEE14-72_contingency.csv"
  genFile = "IEEE14-72_generator.csv"

MyJulia1(rawFile, genFile, contFile)
=#

## Large model
InFile1 = "case.raw"
InFile2 = "case.con"
InFile3 = "case.inl"
InFile4 = "case.rop"

TimeLimitInSeconds = 600
ScoringMethod = 1
NetworkModel = "Network"



include("MyJulia1.jl")
MyJulia1(InFile1, InFile2, InFile3, InFile4, TimeLimitInSeconds, ScoringMethod, NetworkModel)
