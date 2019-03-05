#=
TRASH SCRIPT FOR EVERYTHING
=#

using DataFrames, JuMP

a = Data.generator[:ID][2]

println(a[3:end-1])

println(replace(strip(a), "'" => ""))

b = [parse(Int,replace(strip(x), "'" => "")) for x in Data.generator[:ID]]


a = [1, 2, 3, 1]
b = [4, 5, 6, 4]
c = [(a[i],b[i]) for i = 1:length(a)]


pL = zeros(length(A.caliI),1)
B = Data.load[20 .<= Data.load[:PL] .<= 50, :]
pL[B[:I]] = 100


D_temp = Data.generator[Data.generator[:STAT] .== true, :]

xx = [ v for i in D_temp[:I], v in A.g_temp if i == v[1][1] ]


xx = [parse(Int,replace(strip(x), "'" => "")) for x in Data.generator[:ID]]

yy = colwise(s -> parse(Int,replace(strip(s), "'" => "")), Data.generator[:ID] )

xx = Dict(1 => "ff", 2 => "ddd")


include("Parsing.jl")
using CSV: read

function response_parser!(resp_path, Resp)
  to_skip = 2 # there are to lines at the end "0" and "Q" which indicate the end of the file
  # Resp.data = read(resp_path, footerskip=to_skip, header=names(Resp.data))
  # Resp.data = read(resp_path)
  DD = read(resp_path, DataFrame, delim=',')
  return DD
end

inlFile = "./data/case.inl"
respData = initialize_response_struct()
D = response_parser!(inlFile, respData)

Data = mainData
D = Data.switchedShunt[Data.switchedShunt[:STAT] .== true, :]
rows = size(D)[1]
bCS0 = Dict()
bCS_max = Dict()
bCS_min = Dict()
temp = []
for i = 1:rows
    bCS0[ D[:I][i] ] = D[:BINIT][i] ./ 100
    BL = []
    for j = 1:8
        n = Symbol("N$j")
        b = Symbol("B$j")
        push!(BL, D[n][i] .* D[b][i] ./ 100)
    end
    temp = BL
    NBL = findfirst(BL, 0) - 1
    bCS_max[ D[:I][i] ] = sum([max(0,x) for x in BL])
    bCS_min[ D[:I][i] ] = sum([min(0,x) for x in BL])
end


costsData
rename!(costsData.genDispatch, Dict(:DSPTBL => :TBL))
rename!(costsData.linearCurve, Dict(:LTBL => :CTBL
rename!(costsData.)

D = join(costsData.genDispatch, costsData.activeDispatch, on=[:TBL], kind=:left)

include("components_struct.jl")

d = [1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
bb = Bus(d)

# ==============================

using JuMP, Ipopt
mutable struct G
    id::Tuple{Int64, Int64}
    N::Int64
    G() = new()
end

LG = []
for i=1:2:10, j=1:2
    A = G()
    A.id = (i,j)
    A.N = rand(1:5)
    push!(LG, A)
end

println(LG)
caliG = [g.id for g in LG]

# for i = 1:length(LG)
#     println(LG[i].N)
# end

m = Model(solver=IpoptSolver())

@variable(m, t1[g=1:length(LG), i=1:LG[g].N])



@variable(m, t2[g=LG])

# for g = 1:length(LG)
#     println("g=$g: ", 1:LG[g].N)
# end


for g = 1:length(LG)
    @variable(m, t3[1:LG[g].N])
end

for g in LG
    @variable(m, t4[g.id, 1:g.N])
end

s = 0
for g in LG
    s += g.N
end
print(s)


# ========================
m2 = Model(solver=IpoptSolver())
A = [3,4]
# @macroexpand @variable(m2, x[i=1:length(A), j=1:A[i]])

@variable(m2, y)
x = Array{Array{JuMP.Variable, 1}, 1}()
# @constraintref x_constr =
for i = 1:length(A)
    push!(x, @variable(m2, [1:A[i]], basename="x_$i"))
    # push!(x_constr, @constraint(m2, x[i][1] + 2 == 25))
end
p = PN.GeneratorList[1].p
@variable(m2, p_g[1:6])
@constraint(m2, p' * p_g == 1)

println(m2)




##############################
for g in PN.G
    println(g)
end

#########################################
using DelimitedFiles

path = "./data/case.raw"
rawData = readdlm(path,',', skipblanks=false)

n,m = size(rawData);
println(n, "|", m)

for i = 1:n

    if occursin("END OF LOAD DATA", string(rawData[i,1]),)
        println(rawData[i,1])
        println(i)
    end
end

Data = costsData

D0 = join(Data.genDispatch, Data.activeDispatch, on=[:TBL], kind=:left)
D = join(D0, Data.linearCurve, on=[:CTBL], kind=:left)


# =====================================
using JuMP
using Ipopt

m = Model()
x = @variable(m)
typeof(x)
f = Array{JuMP.VariableRef, 1}()
push!(f, @variable(m, base_name="f_1"))

@constraint(m, con, f[1] >= 0)
