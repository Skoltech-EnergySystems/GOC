include("NetworkData.jl")


#=

=#


"""
constructor for PNetwork

"""
function network_init(PN, Data)
    # SBASE
    PN.s = Data.sbase

    # bus
    PN.caliI = unique(Data.bus[:I])
    PN.caliA = unique(Data.bus[:AREA])
    PN.a = Data.bus[:AREA]
    PN.v0 = Data.bus[:VM]
    PN.theta0 = Data.bus[:VA]
    PN.v_max =Data.bus[:NVHI]
    PN.v_min = Data.bus[:NVLO]
    PN.vK_max = Data.bus[:EVHI]
    PN.vK_min = Data.bus[:EVLO]

    # load
    PN.pL = zeros(length(PN.caliI),1)
    PN.qL = zeros(length(PN.caliI),1)
    # pL qL initialization under condition STATUS == true
    D_temp = Data.load[Data.load[:STATUS] .== true, :]
    PN.pL[D_temp[:I]] += Data.load[:PL] ./ PN.s
    PN.qL[D_temp[:I]] += Data.load[:QL] ./ PN.s

    # fixed shunt data
    # conditional initialization
    PN.gFS = zeros(length(PN.caliI),1)
    PN.bFS = zeros(length(PN.caliI),1)
    # initialization under condition STATUS  == true
    D_temp = Data.fixedBusShunt[Data.fixedBusShunt[:STATUS] .== true, :]
    PN.gFS[D_temp[:I]] += Data.fixedBusShunt[:GL] ./ PN.s
    PN.bFS[D_temp[:I]] += Data.fixedBusShunt[:BL] ./ PN.s

    # generator
    PN.i_g = Data.generator[:I]
    # apply function parse to the column insetead of list comprehensive
    PN.id_g = [parse(Int,replace(strip(x), "'" => "")) for x in Data.generator[:ID]]
    g_temp = [(PN.i_g[j], PN.id_g[j]) for j = 1:length(PN.i_g)]
    PN.caliG = [x for x in Set(g_temp)]
    PN.p0 = Data.generator[:PG] / PN.s
    PN.q0 = Data.generator[:QG] / PN.s
    PN.q_max = Data.generator[:QT] / PN.s
    PN.q_min = Data.generator[:QB] / PN.s
    PN.p_max = Data.generator[:PT] / PN.s
    PN.p_min = Data.generator[:PB] / PN.s

    # init under STAT == true
    D_temp = Data.generator[Data.generator[:STAT] .== true, :]
    PN.G = unique([v for i in D_temp[:I], v in g_temp if i == v[1][1]])


    # [parse(Int,replace(strip(x), "'" => "")) for x in Data.generator[:ID]]
    # line data BRANCH
    PN.iO_e = Data.branch[:I]
    PN.iD_e = Data.branch[:J]
    id_e = [parse(Int,replace(strip(x), "'" => "")) for x in Data.branch[:CKT]]
    e_temp = [(PN.iO_e[j], PN.iD_e[j], id_e[j]) for j = 1:length(PN.iO_e)]
    PN.Eps = unique(e_temp)
    PN.g = Data.branch[:R] ./ (Data.branch[:R] .^ 2 +  Data.branch[:X] .^ 2)
    PN.b = -Data.branch[:X] ./ (Data.branch[:R] .^ 2 + Data.branch[:X] .^ 2)
    PN.bCH = Data.branch[:B]
    PN.R_max = Data.branch[:RATEA] ./ PN.s
    PN.RK_max = Data.branch[:RATEC] ./ PN.s
    # init if ST == 1
    # E
    D_temp = Data.branch[Data.branch[:ST] .== true, :]
    PN.E = unique([v for i in D_temp[:I], v in e_temp if i == v[1][1]])

    # transformer
    PN.iO_f = Data.transformer[:I]
    PN.iD_f = Data.transformer[:J]
    id_f = [parse(Int,replace(strip(x), "'" => "")) for x in Data.transformer[:CKT]]
    f_temp = [(PN.iO_f[j], PN.iD_f[j], id_f[j]) for j = 1:length(PN.iO_f)]
    PN.caliF = unique(f_temp)
    PN.gM_f = Data.transformer[:MAG1]
    PN.bM_f = Data.transformer[:MAG2]
    PN.g_f = Data.transformer[:R12] ./ (Data.transformer[:R12] .^ 2 + Data.transformer[:X12] .^ 2)
    PN.b_f = -Data.transformer[:X12] ./ (Data.transformer[:R12] .^ 2 + Data.transformer[:X12] .^ 2)
    PN.tau_f = Data.transformer[:WINDV1] ./ Data.transformer[:WINDV2]
    PN.theta_f = Data.transformer[:ANG1] .* pi ./ 180
    PN.s_f_max = Data.transformer[:RATA1] ./ PN.s
    PN.sK_f_max = Data.transformer[:RATC1] ./ PN.s

    # conditional initialization
    # F
    D_temp = Data.transformer[Data.transformer[:STAT] .== true, :]
    PN.F = unique([v for i in D_temp[:I], v in f_temp if i == v[1][1]])

    # switched shunt
    # also should be conditional
    """
    parsed row by row of given DataFrame
    should be parsed at the moment of data reading
    """
    D_temp = Data.switchedShunt[Data.switchedShunt[:STAT] .== true, :]
    PN.bCS0, PN.bCS_max, PN.bCS_min  = temp_switchedShunt_init(D_temp, PN.s)

end

function costs_init(Data::CostsStruct, PC::PCosts, s)
    # rename!(Data.genDispatch, Dict(:DSPTBL => :TBL))
    # rename!(Data.linearCurve, Dict(:LTBL => :CTBL))

    D0 = join(Data.genDispatch, Data.activeDispatch, on=[:TBL], kind=:left)
    D = join(D0, Data.linearCurve, on=[:CTBL], kind=:left)

    i = D[:BUS]
    id = D[:GENID]
    # g_temp = [(i[j], id[j]) for j = 1:length(i)]

    caliH = Dict()
    p = Dict()
    c = Dict()
    L = size(D)[1]
    g_ind = []
    for j = 1:L
        g = (D[:BUS][j], D[:GENID][j])
        push!(g_ind, g)
        caliH[g] = D[:NPAIRS][j]
        p[g] = D[:XYmatrix][j][:,1] / s
        c[g] = D[:XYmatrix][j][:,2] / s
    end
    PC.p_gh = p
    PC.c_gh = c
    PC.caliH = caliH
end
