#=
This file containts a structure of The Network with all parameters of the network needed to create a JuMP optimization model

"""
General rules for new variables

    1. \undersocre{x}_i^j translates to xJ_i_min
    2. superscipt indeces become capital letters right afret a variable name
    3. subscript indeces become lower case letters after underscore
    4. all additions like overscore and underscore translate to _min and _max at the end of a      variable name
    5. all caligraphic variables go to cali*, i.e. caliA
"""
=#

# """
# NOT IN USE CURRENTLY
# function to initialize switchedShunt data since it is easier to do line by line
# """
function temp_switchedShunt_init(Data, s)
    BL_len = 8 # the number of (Ni, Bi) pairs
    rows = size(Data)[1]
    bCS0 = Dict()
    bCS_max = Dict()
    bCS_min = Dict()
    for i = 1:rows
        bCS0[ Data[:I][i] ] = Data[:BINIT][i] ./ s
        BL = []
        for j = 1:BL_len
            n = Symbol("N$j")
            b = Symbol("B$j")
            push!(BL, Data[n][i] .* Data[b][i] ./ s)
        end
        NBL = findfirst(BL, 0) - 1
        bCS_max[ Data[:I][i] ] = sum([max(0,x) for x in BL])
        bCS_min[ Data[:I][i] ] = sum([min(0,x) for x in BL])
    end
    return bCS0, bCS_max, bCS_min
end

##############################################
# BUS
mutable struct Bus
    i::Int # i
    a::Int # a
    v0::Float64 # VM
    theta0::Float64 # VA
    v_max::Float64 # NVHI
    v_min::Float64 # NVLO
    vK_max::Float64 # EVHI
    vK_min::Float64 # EVLO

    Bus(Data) = bus_constr(new(), Data)
end
"""
Bus constructor.
in: Data -- line of data from raw file
out: B -- Bus object
"""
function bus_constr(B, Data)     
    B.i = Data[1]
    B.a = Data[5]
    B.v0 = Data[8]
    B.theta0 = Data[9]
    B.v_max = Data[10]
    B.v_min = Data[11]
    B.vK_max = Data[12]
    B.vK_min = Data[13]
    return B
end

"""
Initialize list of Buses from Data
in: Data -- raw data with bus data
    PN -- power network
"""
function Bus_init!(Data, PN)
  n = size(Data)[1]
  PN.BusList = []
  for i = 1:n
    B = Bus(Data[i,:]) # create new bus instance based on current line
    push!(PN.BusList, B) # push new Bus into array of buses
    push!(PN.caliI, B.i) # add its index to set of indeces
    push!(PN.caliA, B.a) # add its area to set of areas
  end
end

# LOAD
mutable struct Load
    i::Int64
    id::Int64
    status::Bool
    pL::Float64
    qL::Float64

    Load(Data, s) = load_constr(new(), Data, s)
end
"""
Load constructor from Data
"""
function load_constr(Ld, Data, s)
    Ld.i = Data[1]
    Ld.id = parse(Int,replace(strip(Data[2]), "'" => ""))
    Ld.status = Data[3]
    Ld.pL = 0.0
    Ld.qL = 0.0
    if Data[3] == true
        Ld.pL += Data[6] ./ s
        Ld.qL += Data[7] ./ s
    end
    return Ld
end
"""
Stores all loads into an array
"""
function Load_init!(Data, PN)
  n = size(Data)[1]
  PN.LoadList = []
  for i = 1:n
    Ld = Load(Data[i,:], PN.s) # new load instance
    push!(PN.LoadList, Ld) # add it to list
  end
end


# FIXED SHUNT
mutable struct FixedShunt
    i::Int
    gFS::Float64
    bFS::Float64
    stat::Bool

    FixedShunt(Data) = fShunt_constr(new(), Data, s)
end
"""
Fixed shunt constructor
"""
function fShunt_constr(fS, Data, s)
    fS.i = Data[1]
    fS.gFS = 0.0
    fS.bFS = 0.0
    fS.stat = Data[3]
    if fS.stat == true
        fS.gFS += Data[4] ./ s
        fS.bFS += Data[5] ./ s
    end
    return fS
end

"""
in challenge 1 it is empty so an empty list created otherwise error "access to undefined reference" occurs
"""
function fShunt_init!(Data, PN)
  n = size(Data)[1]
  PN.fShuntList = []
  for i = 1:n
    B = FixedShunt(Data[i,:], PN.s) # new instance
    push!(PN.fShuntList, B) # stored in list
  end
end

# GENERATOR
mutable struct Generator
    i::Int
    id::Int
    g::Tuple
    p0::Float64
    q0::Float64
    q_max::Float64
    q_min::Float64
    p_max::Float64
    p_min::Float64
    stat::Bool
    # costs
    p::Array
    c::Array
    N::Int
    # participation factor
    alpha::Float64

    Generator(Data, s) = generator_constr(new(), Data, s)
end
"""
Generator constructor
"""
function generator_constr(G, Data, s)
    G.i = Data[1]
    G.id = parse(Int,replace(strip(Data[2]), "'" => ""))
    G.g = (G.i, G.id)
    G.p0 = Data[3] ./ s
    G.q0 = Data[4] ./ s
    G.q_max = Data[5] ./ s
    G.q_min = Data[6] ./ s
    G.p_max = Data[17] ./ s
    G.p_min = Data[18] ./ s
    G.stat = Data[15]
    return G
end

"""
Stores generators into a list of all generators
"""
function Generator_init!(Data, PN)
  n = size(Data)[1]
  PN.GeneratorList = []
  PN.gen_ind_I = Dict{Int, Int}()
  for i = 1:n
    G = Generator(Data[i,:], PN.s) # new generator instance
    push!(PN.GeneratorList, G) # stored into list
    push!(PN.caliG, G.g) # save its compound index in set of generator indeces
    if G.stat == true
        push!(PN.G, G.g) # save its index in set of ACTIVE generator indeces
    end
    # corresondece between compound index g and index i in list of generators for THE SAME GENERATOR
    # generator g can be reached by index PN.gen_ind_I[g] in list of generators PN.GeneratorList
    PN.gen_ind_I[G.g] = length(PN.GeneratorList)
  end
end

# LINE
mutable struct Line
    iO::Int
    iD::Int
    id::Int
    e::Tuple
    g::Float64
    b::Float64
    bCH::Float64
    R_max::Float64
    RK_max::Float64
    st::Bool
    Line(Data, s) = line_constr(new(), Data, s)
end
"""
Line constructor
"""
function line_constr(L, Data, s)
    L.iO = Data[1]
    L.iD = Data[2]
    L.id = parse(Int,replace(strip(Data[3]), "'" => ""))
    L.e = (L.iO, L.iD, L.id)
    L.g = Data[4] / (Data[4] ^ 2 + Data[5] ^ 2)
    L.b = Data[5] / (Data[4] ^ 2 + Data[5] ^ 2)
    L.bCH = Data[6]
    L.R_max = Data[7] / s
    L.RK_max = Data[9] / s
    L.st = Data[14]

    return L
end

"""
Stores lines into an array
"""
function Line_init!(Data, PN)
  n = size(Data)[1]
  PN.LineList = []
  for i = 1:n
    L = Line(Data[i,:], PN.s) # new line instance
    push!(PN.LineList, L) # store into array
    push!(PN.Eps, L.e) # add its index to set of line indeces
    if L.st == true
        push!(PN.E, L.e) # list of active lines
    end
  end
end

# TRANSFORMER
mutable struct Transformer
    iO::Int
    iD::Int
    id::Int
    f::Tuple
    gM::Float64
    bM::Float64
    g::Float64
    b::Float64
    tau::Float64
    theta::Float64
    s_max::Float64
    sK_max::Float64
    st::Bool

    Transformer(Data, s) = transformer_constr(new(), Data, s)
end
"""
Transformer constructor
"""
function transformer_constr(T, Data, s)
    T.iO = Data[1]
    T.iD = Data[2]
    T.id = parse(Int,replace(strip(Data[4]), "'" => ""))
    T.f = (T.iO, T.iD, T.id)
    T.gM = Data[8]
    T.bM = Data[9]
    T.g = Data[22] / (Data[22] ^ 2 + Data[23] ^ 2)
    T.b = Data[23] / (Data[22] ^ 2 + Data[23] ^ 2)
    T.tau = Data[25] ./ Data[42]
    T.theta = Data[27] .* pi ./ 180
    T.s_max = Data[28] ./ s
    T.sK_max = Data[30] ./ s
    T.st = Data[12]

    return T
end

"""
Stores transformets into array
"""
function Transformer_init!(Data, PN)
  n = size(Data)[1]
  PN.TransformerList = []
  step = 4 # these 4 lines should be contacatinated as one to put into the frame
  lines_len = [21, 3, 17, 2] # length of each of 4 lines. 43 attributes of transformer in total
  for i = 1:step:n
    full_line = []
    for j = 0:step-1
      full_line = [full_line; Data[i+j, 1:lines_len[j+1]]]; # combined line
    end
    T = Transformer(full_line, PN.s) # new Transfromer instance
    push!(PN.TransformerList, T) # store into array
    push!(PN.caliF, T.f) # save its index
    if T.st == true
        push!(PN.F, T.f) # save to set of active
    end
  end
end

# SWITCHED SHUNT
mutable struct SwitchedShunt
    i::Int
    bCS0::Float64
    bCS_max::Float64
    bCS_min::Float64
    stat::Bool

    SwitchedShunt(Data, s) = sShunt_constr(new(), Data, s)
end
"""
SwitchedShunt constructor
"""
function sShunt_constr(sS, Data, s)
    sS.i = Data[1]
    BL_len = 8
    BL = []
    st_ind = 11 # index of N1, starting the sequence of Ni, Bi
    # all that is according to the instructions in the GOC file
    if Data[4] == true # if active
        sS.stat = true
        sS.bCS0 = Data[10]
        for j = st_ind:2:2*BL_len-1
            push!(BL, Data[j] * Data[j+1] / s)
        end
        NBL = findfirst(isequal(0), BL) - 1
        sS.bCS_max = sum([max(0,x) for x in BL])
        sS.bCS_min = sum([min(0,x) for x in BL])
    end
    return sS
end
"""
Stores SwitchedShunts into array
"""
function sShunt_init!(Data, PN)
  n = size(Data)[1]
  PN.sShuntList = []
  PN.get_shunt_index = Dict{Int, Int}()
  for i = 1:n
    S = SwitchedShunt(Data[i,:], PN.s)
    push!(PN.sShuntList, S)
    PN.get_shunt_index[S.i] = length(PN.sShuntList)
  end
end

#########################################
# Power Network structure
mutable struct PNetwork
    # base MVA
    s::Float16

    #bus data
    caliI::Set # bus indeces
    caliA::Set # areas
    BusList::Array{Bus, 1} # array of Buses

    # load data
    LoadList::Array{Load, 1} # array of Loads

    # fixed shunt data
    fShuntList::Array{FixedShunt}  # array of Fixed Shunts

    # generator data
    caliG::Set # set of generators
    G::Set # set of active generators
    GeneratorList::Array{Generator} # array of Generators
    # correspondance between index g and and index in the list
    gen_ind_I::Dict{Tuple, Int} # change to get_gener_index

    # line data BRANCH
    Eps::Set # lines' indeces
    E::Set # active lines' indeces
    LineList::Array{Line} # array of Lines

    # transformer data
    caliF::Set # transformers' indeces
    F::Set # active transformers' indeces
    TransformerList::Array{Transformer} # array of Transformers

    # switched shunt
    sShuntList::Array{SwitchedShunt} # array of Switched Shunts
    # corresondance between switched shunt index e and its index in the list
    get_shunt_index::Dict{Int, Int}

    PNetwork() = constr_network(new())
end
"""
Power Network constructor
"""
function constr_network(PN)
    PN.caliI = Set(Array{Int64, 1}())
    PN.caliA = Set(Array{Int64, 1}())
    PN.caliG = Set(Array{Tuple{Int64, Int64}, 1}())
    PN.G = Set(Array{Tuple{Int64, Int64}, 1}())
    PN.Eps = Set(Array{Tuple{Int64, Int64, Int64}, 1}())
    PN.E = Set(Array{Tuple{Int64, Int64, Int64}, 1}())
    PN.caliF = Set(Array{Tuple{Int64, Int64, Int64}, 1}())
    PN.F = Set(Array{Tuple{Int64, Int64, Int64}, 1}())

    return PN
end

# =====================================
# Costs data
"""
Initializes costs data for Power Network PN
in: Data -- DataFrames with all costs
    PN -- Power Network
"""
function costs_init!(Data, PN)
    # join DataFrames as described in the instructions file
    D0 = join(Data.genDispatch, Data.activeDispatch, on=[:TBL], kind=:left)
    D = join(D0, Data.linearCurve, on=[:CTBL], kind=:left)
    # parse line by line
    for i = 1:size(D)[1]
        L = D[i,:] # next line
        g = (L[:BUS][1], L[:GENID][1]) # generator's index g
        j = PN.gen_ind_I[g] # its index in the list

        PN.GeneratorList[j].N = L[:NPAIRS][1] # number of segments of piecewise linear costs
        PN.GeneratorList[j].p = L[:XYmatrix][:,1] / PN.s # generation points
        PN.GeneratorList[j].c = L[:XYmatrix][:,2] / PN.s # costs
    end
end

###################################
# OLD STRUCTURE TO DELETE
mutable struct PCosts
    caliH::Dict
    p_gh::Dict
    c_gh::Dict
    PCosts() = new()
end
