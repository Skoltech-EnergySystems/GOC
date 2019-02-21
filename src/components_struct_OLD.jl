"""
copied to NetworkData.jl
"""


# BUS
mutable struct Bus
    # caliI::Int
    # caliA::Int
    a::Int # a
    v0::Float64 # VM
    theta0::Float64 # VA
    v_max::Float64 # NVHI
    v_min::Float64 # NVLO
    vK_max::Float64 # EVHI
    vK_min::Float64 # EVLO

    # Bus() = new()
    Bus(Data) = bus_constr(new(), Data)
end
function bus_constr(B, Data)
    # B = new()
    # B.a = Data.bus[:AREA]
    # B.v0 = Data.bus[:VM]
    # B.theta0 = Data.bus[:VA]
    # B.v_max =Data.bus[:NVHI]
    # B.v_min = Data.bus[:NVLO]
    # B.vK_max = Data.bus[:EVHI]
    # B.vK_min = Data.bus[:EVLO]
    B.a = Data[5]
    B.v0 = Data[8]
    B.theta0 = Data[9]
    B.v_max = Data[10]
    B.v_min = Data[11]
    B.vK_max = Data[12]
    B.vK_min = Data[13]
    return B
end

# LOAD
mutable struct Load
    pL::Float64
    qL::Float64

    Load(Data, s) = load_constr(new(), Data, s)
end
function load_constr(Ld, Data, s)
    Ld.pL = 0.0
    Ld.qL = 0.0
    if Data[3] == true
        Ld.pL += Data[6] ./ s
        Ld.qL += Data[7] ./ s
    end
    return Ld
end

# FIXED SHUNT
mutable struct FixedShunt
    gFS::Float64
    bFS::Float64

    FixedShunt(Data) = fShunt_constr(new(), Data, s)
end
function fShunt_constr(fS, Data, s)
    fS.gFS = 0.0
    fS.bFS = 0.0
    if Data[3] == true
        fS.gFS += Data[4] ./ s
        fS.bFS += Data[5] ./ s
    end
    return fS
end

# GENERATOR
mutable struct Generator
    # caliG::Int
    # G::Int
    i::Int
    id::Int
    g::Tuple
    p0::Float64
    q0::Float64
    q_max::Float64
    q_min::Float64
    p_max::Float64
    p_min::Float64
    # costs
    p::Array
    c::Array
    N::Int
    alpha::Float64 # participation factor

    Generator(Data, s) = generator_constr(new(), Data, s)
end
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
    return G
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

    Line(Data, s) = line_constr(new(), Data, s)
end
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

    return L
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

    Transformer(Data, s) = transformer_constr(new(), Data, s)
end
function transformer_constr(T, Data, s)
    T.iO = Data[1]
    T.iD = Data[2]
    T.id = parse(Int,replace(strip(Data[4]), "'" => ""))
    T.f = (T.iO, T.iD, T.id)
    T.gM = Data[8]
    T.bM = Data[9]
    # println("R12 $(typeof(Data[22])), | X12 $(typeof(Data[23]))")
    T.g = Data[22] / (Data[22] ^ 2 + Data[23] ^ 2)
    T.b = Data[23] / (Data[22] ^ 2 + Data[23] ^ 2)
    T.tau = Data[25] ./ Data[42]
    T.theta = Data[27] .* pi ./ 180
    T.s_max = Data[28] ./ s
    T.sK_max = Data[30] ./ s

    return T
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
function sShunt_constr(sS, Data, s)
    sS.i = Data[1]
    BL_len = 8
    BL = []
    st_ind = 11 # index of N1, starting the sequence of Ni, Bi
    if Data[4] == true
        sS.stat = true
        sS.bCS0 = Data[10]
        for j = st_ind:2:2*BL_len-1
            push!(BL, Data[j] * Data[j+1] / s)
        end
        NBL = findfirst(BL, 0) - 1
        sS.bCS_max = sum([max(0,x) for x in BL])
        sS.bCS_min = sum([min(0,x) for x in BL])
    end
    return sS
end
