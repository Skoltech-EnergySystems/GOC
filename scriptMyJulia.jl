### Script
include("MyJulia1.jl");
# initial model (small)
##=
 rawFile = "IEEE14-1_powersystem.raw";
 contFile = "IEEE14-1_contingency.csv";
 genFile = "IEEE14-1_generator.csv";

MyJulia1(rawFile, genFile, contFile)

#=
## Large model
  rawFile = "RTS96-1_powersystem.raw"
  contFile = "RTS96-1_contingency.csv"
  genFile = "RTS96-1_generator.csv"

  MyJulia1(rawFile, genFile, contFile)
=#
