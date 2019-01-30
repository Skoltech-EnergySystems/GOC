
#=
"""
General rules for new variables

    1. \undersocre{x}_i^j translates to xJ_i_min
    2. superscipt indeces become capital letters right afret a variable name
    3. subscript indeces become lower case letters after underscore
    4. all additions like overscore and underscore translate to _min and _max at the end of a      variable name
    5. all caligraphic variables go to cali*, i.e. caliA
"""
=#

# main data descrbing network parameters
mutable struct PNetwork

    # main sets
    caliH::Array # costs
    caliK::Array # contingencies

    # base MVA
    s::Float16

    #bus data
    # i = bus[:I] a = bus[:AREA]
    caliI::Array
    caliA::Array
    a::Array # a
    v0::Array # VM
    theta0::Array # VA
    v_max::Array # NVHI
    v_min::Array # NVLO
    vK_max::Array # EVHI
    vK_min::Array # EVLO

    # load data
    pL::Array
    qL::Array

    # fixed shunt data
    gFS::Array
    bFS::Array

    # generator data
    # i = generator[:I], id = generator[:ID], g = (i, id)
    caliG::Array
    G::Array
    i_g::Array
    id_g::Array
    p0::Array
    q0::Array
    q_max::Array
    q_min::Array
    p_max::Array
    p_min::Array

    # line data BRANCH
    # i = line[:I] iPrime = line[:J] id = line[:CKT] e = (i, iPrime, id)
    Eps::Array
    E::Array
    iO_e::Array
    iD_e::Array
    g::Array
    b::Array
    bCH::Array
    R_max::Array
    RK_max::Array

    # transformer data
    # i = :I, iPrime = :J, id = :CKT, f = (i, iPrime, id)
    caliF::Array
    F::Array
    iO_f::Array
    iD_f::Array
    gM_f::Array
    bM_f::Array
    g_f::Array
    b_f::Array
    tau_f::Array
    theta_f::Array
    s_f_max::Array
    sK_f_max::Array

    # swtiched shunt
    """
    needs to be reformulated
    """
    bCS0::Dict
    # BL::Array
    # BL_8::Float16
    # NBL::Float16
    bCS_max::Dict
    bCS_min::Dict

    # to be continued...

    """
    TO DO
    write custom constructor as
    PNetwork() = some_function(new(), Data)
    """
    PNetwork() = new()
end



"""
function to initialize switchedShunt data since it is easier to do line by line
"""
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



# =====================================
# Costs data
mutable struct PCosts
    caliH::Dict
    p_gh::Dict
    c_gh::Dict
    PCosts() = new()
end
