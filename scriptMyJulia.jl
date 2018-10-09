### Script
include("MyJulia1.jl");
cd("C:/Users/Ильгиз/Documents/Документы_Ильгиз/Skoltech_2018/Grid Competition/Code/GOC_v14_added_multiple schemes")
# initial model (small)
##=
 rawFile = "IEEE14-1_powersystem.raw";
 contFile = "IEEE14-1_contingency.csv";
 genFile = "IEEE14-1_generator.csv";

MyJulia1(rawFile, genFile, contFile)
#=
# 14 nodes but 72 contingencies
  rawFile = "IEEE14-72_powersystem.raw"
  contFile = "IEEE14-72_contingency.csv"
  genFile = "IEEE14-72_generator.csv"

  MyJulia1(rawFile, genFile, contFile)
=#
#=
## Large model
  rawFile = "RTS96-1_powersystem.raw"
  contFile = "RTS96-1_contingency.csv"
  genFile = "RTS96-1_generator.csv"

  MyJulia1(rawFile, genFile, contFile)
=#
