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



using JuMP
mutable struct G
    g::Set
    N::Int16
    G() = new(g, N)
end

LG = []
for i=1:2:10, j=1:2
    G((i,j), rand(1:5, 1))
    push!(LG, G)
end

println(rand(1:5, 1))
