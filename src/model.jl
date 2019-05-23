using JuMP
using Ipopt
include("NetworkData.jl")
include("DataStructures.jl")

function create_model(PN::PNetwork, continData::ContingenciesStruct)

    # OPF = Model(solver=IpoptSolver())
    OPF = Model(with_optimizer(Ipopt.Optimizer))


    # base modell

    # """
    # loop through all complex indeces like g = (i, id) and add constraints manually.
    # """
    # base
    # @VariableRefs OPF begin
    #     c
    #     c_g
    #     cSigma
    #     t_gh
    #     v_i
    #     theta_i
    #     bCS_i
    #     sigPp_i
    #     sigPm_i
    #     sigQp_i
    #     sigQm_i
    #     sigS_e
    #     sigS_f
    #     p_g
    #     q_g
    #     pO_e
    #     pD_e
    #     qO_e
    #     qD_e
    #     pO_f
    #     pD_f
    #     qO_f
    #     qD_f
    # end

    # contingency
    # @VariableRefs OPF begin
    #     cSig_k
    #     delta_k
    #     v_ik
    #     theta_ik
    #     bCS_ik
    #     sigPp_ik
    #     sigPm_ik
    #     sigQp_ik
    #     sigQm_ik
    #     sigS_ek
    #     sigS_fk
    #     p_gk
    #     q_gk
    #     pO_ek
    #     pD_ek
    #     qO_ek
    #     qD_ek
    #     pO_fk
    #     pD_fk
    #     qO_fk
    #     qD_fk
    # end

# ∀ g ∈ G
    p_g = Array{JuMP.VariableRef, 1}()
    q_g = Array{JuMP.VariableRef, 1}()
    t_gh = Array{Array{JuMP.VariableRef, 1}, 1}()
    c_g = Array{JuMP.VariableRef, 1}()

    lenG = length(PN.G)
    p_g_constr = Array{JuMP.ConstraintRef, 1}() # 33-34
    q_g_constr = Array{JuMP.ConstraintRef, 1}() # 35-36
    c_g_constr = Array{JuMP.ConstraintRef, 1}() # 2
    t_g_constr = Array{JuMP.ConstraintRef, 1}() # 5
    diffG = setdiff(PN.caliG, PN.G)

    sizeG = length(PN.caliG);
    for i in 1:sizeG # index of generator[g] in the list of all generators
# 33
        # push!(p_g, @variable(OPF, lower_bound=PN.GeneratorList[i].p_min, upper_bound=PN.GeneratorList[i].p_max, base_name="p_$i"))
# 35
        # push!(q_g, @variable(OPF, lower_bound=PN.GeneratorList[i].q_min, upper_bound=PN.GeneratorList[i].q_max, base_name="q_$i"))
# 4, 5
        push!(t_gh, @variable(OPF, [h=1:PN.GeneratorList[i].N], lower_bound=0, base_name="t_$i"))
        push!(t_g_constr, @constraint(OPF, sum(t_gh[end]) == 1))
# 34, 36
        if i in diffG
            push!(q_g, @variable(OPF, base_name="q_$i"))
            push!(p_g, @variable(OPF, base_name="p_$i"))
            push!(p_g_constr, @constraint(OPF, p_g[end] == 0))
            push!(q_g_constr, @constraint(OPF, q_g[end] == 0))
        else
            push!(q_g, @variable(OPF, lower_bound=PN.GeneratorList[i].q_min, upper_bound=PN.GeneratorList[i].q_max, base_name="q_$i"))
            push!(p_g, @variable(OPF, lower_bound=PN.GeneratorList[i].p_min, upper_bound=PN.GeneratorList[i].p_max, base_name="p_$i"))
# 3
            push!(p_g_constr, @constraint(OPF, p_g[end] == PN.GeneratorList[i].p' * t_gh[end]))
        end
# 2
        push!(c_g, @variable(OPF, base_name="c_$i"))
        push!(c_g_constr, @constraint(OPF, c_g[end] == PN.GeneratorList[i].c' * t_gh[end]))
    end

# ∀ i ∈ I
    v_i = Array{JuMP.VariableRef, 1}()
    theta_i = Array{JuMP.VariableRef, 1}()
    bCS_i = Array{JuMP.VariableRef, 1}()
    sigPp_i = Array{JuMP.VariableRef, 1}()
    sigPm_i = Array{JuMP.VariableRef, 1}()
    sigQp_i = Array{JuMP.VariableRef, 1}()
    sigQm_i = Array{JuMP.VariableRef, 1}()

    sigP_diff = Array{JuMP.ConstraintRef, 1}() # 46
    sigQ_diff = Array{JuMP.ConstraintRef, 1}() # 49

    for i in PN.caliI
# 32
        push!(v_i, @variable(OPF, lower_bound=PN.BusList[i].v_min, upper_bound=PN.BusList[i].v_max, base_name="v_$i"))

# 37
        if i in keys(PN.get_sShunt_index)
            push!(bCS_i, @variable(OPF, lower_bound=PN.sShuntList[PN.get_sShunt_index[i]].bCS_min, upper_bound=PN.sShuntList[PN.get_sShunt_index[i]].bCS_max, base_name="bCS_$(PN.get_sShunt_index[i])"))
        end
# 47 48 50 51
        push!(sigPp_i, @variable(OPF, lower_bound=0, base_name="sigPp_$i") )
        push!(sigPm_i, @variable(OPF, lower_bound=0, base_name="sigPm_$i") )
        push!(sigQp_i, @variable(OPF, lower_bound=0, base_name="sigQp_$i") )
        push!(sigQm_i, @variable(OPF, lower_bound=0, base_name="sigQm_$i") )
# 46
        # push!(sigP_diff, @constraint(OPF, sigPp_i[end] - sigPm_i[end] == )  )

    end
#
    sizeI = maximum(PN.caliI);
    @variable(OPF, theta_i[1:sizeI])
    @variable(OPF, v_i[1:sizeI])

    # Line flow: 38 - 41
    pO_e = Array{JuMP.VariableRef, 1}()
    qO_e = Array{JuMP.VariableRef, 1}()
    pD_e = Array{JuMP.VariableRef, 1}()
    qD_e = Array{JuMP.VariableRef, 1}()

    pO_e_constr = Array{JuMP.ConstraintRef, 1}()
    qO_e_constr = Array{JuMP.ConstraintRef, 1}()
    pD_e_constr = Array{JuMP.ConstraintRef, 1}()
    qD_e_constr = Array{JuMP.ConstraintRef, 1}()
    # ∀ e ∈ E
    for L in PN.LineList
        push!(pO_e, @variable(OPF, base_name="pO_e_$(L.e)"))
        push!(qO_e, @variable(OPF, base_name="qO_e_$(L.e)"))
        push!(pD_e, @variable(OPF, base_name="pD_e_$(L.e)"))
        push!(qD_e, @variable(OPF, base_name="qD_e_$(L.e)"))
        # 38
        push!(pO_e_constr, @NLconstraint(OPF,
        pO_e[end] == L.g * v_i[L.iO]^2 + (-L.g * cos(theta_i[L.iO] - theta_i[L.iD]) - L.b * sin(theta_i[L.iO] - theta_i[L.iD])) * v_i[L.iO] * v_i[L.iD] ))
        # 39
        push!(qO_e_constr, @NLconstraint(OPF,
        qO_e[end] == -(L.b + L.bCH/2) * v_i[L.iO]^2 + (L.b * cos(theta_i[L.iO] - theta_i[L.iD]) - L.g * sin(theta_i[L.iO] - theta_i[L.iD])) * v_i[L.iO] * v_i[L.iD] ))
        # 40
        push!(pD_e_constr, @NLconstraint(OPF,
        pD_e[end] == L.g * v_i[L.iD]^2 + (-L.g * cos(theta_i[L.iD] - theta_i[L.iO]) - L.b * sin(theta_i[L.iD] - theta_i[L.iO])) * v_i[L.iO] * v_i[L.iD] ))
        # 41
        push!(qD_e_constr, @NLconstraint(OPF,
        qD_e[end] == -(L.b + L.bCH/2) * v_i[L.iD]^2 + (L.b * cos(theta_i[L.iD] - theta_i[L.iO]) - L.g * sin(theta_i[L.iD] - theta_i[L.iO])) * v_i[L.iO] * v_i[L.iD] ))
    end

    # Transformer flow: 42 - 45
    pO_f = Array{JuMP.VariableRef, 1}()
    qO_f = Array{JuMP.VariableRef, 1}()
    pD_f = Array{JuMP.VariableRef, 1}()
    qD_f = Array{JuMP.VariableRef, 1}()

    pO_f_constr = Array{JuMP.ConstraintRef, 1}()
    qO_f_constr = Array{JuMP.ConstraintRef, 1}()
    pD_f_constr = Array{JuMP.ConstraintRef, 1}()
    qD_f_constr = Array{JuMP.ConstraintRef, 1}()
    #
    # ∀ f ∈ F
    for T in PN.TransformerList
        push!(pO_f, @variable(OPF, base_name="pO_f_$(T.f)"))
        push!(qO_f, @variable(OPF, base_name="qO_f_$(T.f)"))
        push!(pD_f, @variable(OPF, base_name="pD_f_$(T.f)"))
        push!(qD_f, @variable(OPF, base_name="qD_f_$(T.f)"))
        #
        # 42
        push!(pO_f_constr, @NLconstraint(OPF,
        pO_f[end] == (T.g / T.tau^2 + T.gM) * v_i[T.iO]^2 + (-T.g / T.tau * cos(theta_i[T.iO] - theta_i[T.iD]  - T.theta) - T.b / T.tau * sin(theta_i[T.iO] - theta_i[T.iD])) * v_i[T.iO] * v_i[T.iD] ))
        # 43
        push!(qO_f_constr, @NLconstraint(OPF,
        qO_f[end] == -(T.b / T.tau^2 + T.bM) * v_i[T.iO]^2 + (T.b / T.tau * cos(theta_i[T.iO] - theta_i[T.iD]) - T.g / T.tau * sin(theta_i[T.iO] - theta_i[T.iD])) * v_i[T.iO] * v_i[T.iD] ))
        # 44
        push!(pD_f_constr, @NLconstraint(OPF,
        pD_f[end] == T.g * v_i[T.iD]^2 + (-T.g / T.tau * cos(theta_i[T.iD] - theta_i[T.iO]) - T.b / T.tau * sin(theta_i[T.iD] - theta_i[T.iO])) * v_i[T.iO] * v_i[T.iD] ))
        # 45
        push!(qD_f_constr, @NLconstraint(OPF,
        qD_f[end] == -T.b * v_i[T.iD]^2 + (T.b / T.tau * cos(theta_i[T.iD] - theta_i[T.iO]) - T.g / T.tau * sin(theta_i[T.iD] - theta_i[T.iO])) * v_i[T.iO] * v_i[T.iD] ))
    end

    # Nodal equations
    # 46
    pNod_constr = Array{JuMP.ConstraintRef, 1}()
    for i = 1:sizeI
        push!(pNod_constr, @NLconstraint(OPF,
            sum(p_g[g] for g in 1:sizeG if PN.GeneratorList[g].i==i) - sum(PN.LoadList[l].pL for l in 1:length(PN.LoadList) if PN.LoadList[l].i==i)
            - sum(PN.fShuntList[f].gFS*v_i[f]^2 for f in 1:length(PN.fShuntList) if PN.fShuntList[f].i==i) # no fixed shunts in Challenge 1
            - sum(pO_e[oe] for oe in 1:length(PN.LineList) if PN.LineList[oe].iO==i) - sum(pD_e[de] for de in 1:length(PN.LineList) if PN.LineList[de].iD==i)
            - sum(pO_f[of] for of in 1:length(PN.TransformerList) if PN.TransformerList[of].iO==i)
            - sum(pD_f[df] for df in 1:length(PN.TransformerList) if PN.TransformerList[df].iD==i)
            == sigPp_i[i] - sigPm_i[i]))
    end

    # 49
    qNod_constr = Array{JuMP.ConstraintRef, 1}()
    for i = 1:sizeI
        push!(qNod_constr, @NLconstraint(OPF,
            sum(q_g[g] for g in 1:sizeG if PN.GeneratorList[g].i==i) - sum(PN.LoadList[l].qL for l in 1:length(PN.LoadList) if PN.LoadList[l].i==i)
            - sum(-PN.sShuntList[s].bCS0*v_i[s]^2 for s in 1:length(PN.sShuntList) if PN.sShuntList[s].i==i) # no fixed shunts in Challenge 1
            - sum(qO_e[oe] for oe in 1:length(PN.LineList) if PN.LineList[oe].iO==i) - sum(qD_e[de] for de in 1:length(PN.LineList) if PN.LineList[de].iD==i)
            - sum(qO_f[of] for of in 1:length(PN.TransformerList) if PN.TransformerList[of].iO==i)
            - sum(qD_f[df] for df in 1:length(PN.TransformerList) if PN.TransformerList[df].iD==i)
            == sigQp_i[i] - sigQm_i[i]))
    end

    # ∀ e ∈ E
    # Line current ratings: 52-54
    sigS_e = Array{JuMP.VariableRef, 1}()
    CurrentO_e_constr = Array{JuMP.ConstraintRef, 1}()
    CurrentD_e_constr = Array{JuMP.ConstraintRef, 1}()
    for L in PN.LineList
        # 53
        push!(sigS_e, @variable(OPF, lower_bound=0, base_name="sigS_$(L.e)") )
        # 52
        push!(CurrentO_e_constr, @NLconstraint(OPF, sqrt(pO_e[end]^2 + qO_e[end]^2) ≤ L.R_max*v_i[L.iO] + sigS_e[end]))
        # 54
        push!(CurrentD_e_constr, @NLconstraint(OPF, sqrt(pD_e[end]^2 + qD_e[end]^2) ≤ L.R_max*v_i[L.iD] + sigS_e[end]))
    end

    # ∀ f ∈ F
    # Transformer power ratings: 55-57
    sigS_f = Array{JuMP.VariableRef, 1}()
    PowerO_f_constr = Array{JuMP.ConstraintRef, 1}()
    PowerD_f_constr = Array{JuMP.ConstraintRef, 1}()
    for T in PN.TransformerList
        # 56
        push!(sigS_f, @variable(OPF, lower_bound=0, base_name="sigS_$(T.f)") )
        # 55
        push!(PowerO_f_constr, @NLconstraint(OPF, sqrt(pO_f[end]^2 + qO_f[end]^2) <= T.s_max + sigS_f[end]))
        # 57
        push!(PowerD_f_constr, @NLconstraint(OPF, sqrt(pD_f[end]^2 + qD_f[end]^2) <= T.s_max + sigS_f[end]))
    end

    ####CONTINGENCIES
    ###GENERATORS CONTINGENCIES

    p_g_kg = Array{JuMP.VariableRef, 1}()
    q_g_kg = Array{JuMP.VariableRef, 1}()
    # t_gh_kg = Array{Array{JuMP.VariableRef, 1}, 1}()
    c_g_kg = Array{JuMP.VariableRef, 1}()

    # lenG = length(PN.G)
    p_g_constr_kg = Array{JuMP.ConstraintRef, 1}() # 33-34
    q_g_constr_kg = Array{JuMP.ConstraintRef, 1}() # 35-36
    c_g_constr_kg = Array{JuMP.ConstraintRef, 1}() # 2
    # t_g_constr_kg = Array{JuMP.ConstraintRef, 1}() # 5


    v_ik = Array{JuMP.VariableRef, 1}()
    theta_ik = Array{JuMP.VariableRef, 1}()
    bCS_i_kg = Array{JuMP.VariableRef, 1}()
    sigPp_i_kg = Array{JuMP.VariableRef, 1}()
    sigPm_i_kg = Array{JuMP.VariableRef, 1}()
    sigQp_i_kg = Array{JuMP.VariableRef, 1}()
    sigQm_i_kg = Array{JuMP.VariableRef, 1}()


    sigP_diff_kg = Array{JuMP.ConstraintRef, 1}() # 46
    sigQ_diff_kg = Array{JuMP.ConstraintRef, 1}() # 49

    caliK = length(continData.genCont)
    for j in 1:caliK
        bad_generator = (parse(Int64, continData.genCont[j]), 1)
        lenG = length(PN.G)
        PN.G  = setdiff(PN.G, Set([bad_generator])) # Generators that are online in the contringency case
        diffG = setdiff(PN.caliG, PN.G) # Generators that are switched in the contingency case
        # diffG = push!(diffG, bad_generator)


        sizeG = length(PN.caliG);
        for i in 1:sizeG # index of generator[g] in the list of all generators

    # 4, 5
            # push!(t_gh_kg, @variable(OPF, [h=1:PN.GeneratorList[i].N], lower_bound=0, base_name="t_$i"*"_$j"*"_g"))
            # push!(t_g_constr_kg, @constraint(OPF, sum(t_gh_kg[end]) == 1))
    # 60, 62

    ##DONE Above!!!!!!!!!!!!!!
            if i in diffG
                push!(p_g_kg, @variable(OPF, base_name="q_$i"*"_$j"*"_g"))
                push!(q_g_kg, @variable(OPF, base_name="q_$i"*"_$j"*"_g"))
                push!(p_g_constr_kg, @constraint(OPF, p_g_kg[end] == 0))
                push!(q_g_constr_kg, @constraint(OPF, q_g_kg[end] == 0))
            else
        # 59
                push!(p_g_kg, @variable(OPF, lower_bound=PN.GeneratorList[i].p_min, upper_bound=PN.GeneratorList[i].p_max, base_name="q_$i"*"_$j"*"_g"))
        # 61
                push!(q_g_kg, @variable(OPF, lower_bound=PN.GeneratorList[i].q_min, upper_bound=PN.GeneratorList[i].q_max, base_name="q_$i"*"_$j"*"_g"))
    # 3
                # push!(p_g_constr_kg, @constraint(OPF, p_g_kg[end] == PN.GeneratorList[i].p' * t_gh_kg[end]))
            end
    # DONE
            # push!(c_g_kg, @variable(OPF, base_name="c_$i"*"_$j"*"_g"))
            # push!(c_g_constr_kg, @constraint(OPF, c_g_kg[end] == PN.GeneratorList[i].c' * t_gh_kg[end]))
        end

    # ∀ i ∈ I, j ∈ J

    #DONE !!!!!!! ABOVE
        for i in PN.caliI
    # 58
            push!(v_ik, @variable(OPF, lower_bound=PN.BusList[i].vK_min, upper_bound=PN.BusList[i].vK_max, base_name="v_$i"*"_$j"*"_g"))
            push!(theta_ik, @variable(OPF, base_name="theta_$i"*"_$j"*"_g"))

    # 63
            if i in keys(PN.get_sShunt_index)
                # push!(bCS_i_kg, @variable(OPF, lower_bound=PN.sShuntList[PN.get_sShunt_index[i]].bCS_min, upper_bound=PN.sShuntList[PN.get_sShunt_index[i]].bCS_max, base_name="bCS_$(PN.get_sShunt_index[i])"*"_$j_g"))
                push!(bCS_i_kg, @variable(OPF, lower_bound=PN.sShuntList[PN.get_sShunt_index[i]].bCS_min, upper_bound=PN.sShuntList[PN.get_sShunt_index[i]].bCS_max))
            end
    # 73 74 76 77
            #DONE ABOVE
            push!(sigPp_i_kg, @variable(OPF, lower_bound=0, base_name="sigPp_$i"*"_$j"*"_g") )
            push!(sigPm_i_kg, @variable(OPF, lower_bound=0, base_name="sigPm_$i"*"_$j"*"_g") )
            push!(sigQp_i_kg, @variable(OPF, lower_bound=0, base_name="sigQp_$i"*"_$j"*"_g") )
            push!(sigQm_i_kg, @variable(OPF, lower_bound=0, base_name="sigQm_$i"*"_$j"*"_g") )
    #
            # push!(sigP_diff, @constraint(OPF, sigPp_i[end] - sigPm_i[end] == )  )

        end

        #DONE ABOVE. Waintig for 46 and 49 constraints.
        sizeI = maximum(PN.caliI);
        # @variable(OPF, theta_i[1:sizeI])
        # @variable(OPF, v_i[1:sizeI])
    end

    # Line flow: 78 - 80
    pO_e_kl = Array{JuMP.VariableRef, 1}()
    qO_e_kl = Array{JuMP.VariableRef, 1}()
    pD_e_kl = Array{JuMP.VariableRef, 1}()
    qD_e_kl = Array{JuMP.VariableRef, 1}()

    pO_e_constr_kl = Array{JuMP.ConstraintRef, 1}()
    qO_e_constr_kl = Array{JuMP.ConstraintRef, 1}()
    pD_e_constr_kl = Array{JuMP.ConstraintRef, 1}()
    qD_e_constr_kl = Array{JuMP.ConstraintRef, 1}()


    #for j in 1:length(continData.genCont)
    #end


    # ∀ e ∈ E
    # Line current ratings:78-80
    LineList_cont = [];
    for k in 1:length(continData.lineCont)

        LineList_k = PN.LineList
        left = continData.lineCont[k][1]
        right = continData.lineCont[k][2]
        for k in 1:length(PN.LineList)
            if PN.LineList[k].iO == left && PN.LineList[k].iD == right
                LineList_cont = [LineList_cont;PN.LineList[k:k]]
                # LineList_k = setdiff(PN.LineList, PN.LineList[k:k])
            end
        end
        # print(size(LineList_cont)) # Size is correct

        LineList_k = setdiff(PN.LineList, LineList_cont); # Lines that are active in the contingency case

        for L in LineList_k
            # print(size(LineList_k)) # Size is correct
            push!(pO_e_kl, @variable(OPF, base_name="pO_e_$(L.e)_kl"))
            push!(qO_e_kl, @variable(OPF, base_name="qO_e_$(L.e)_kl"))
            push!(pD_e_kl, @variable(OPF, base_name="pD_e_$(L.e)_kl"))
            push!(qD_e_kl, @variable(OPF, base_name="qD_e_$(L.e)_kl"))

            # 64
            push!(pO_e_constr_kl, @NLconstraint(OPF,
            pO_e_kl[end] == L.g * v_ik[L.iO]^2 + (-L.g * cos(theta_ik[L.iO] - theta_ik[L.iD]) - L.b * sin(theta_ik[L.iO] - theta_ik[L.iD])) * v_ik[L.iO] * v_ik[L.iD] ))
            # 65
            push!(qO_e_constr_kl, @NLconstraint(OPF,
            qO_e_kl[end] == -(L.b + L.bCH/2) * v_ik[L.iO]^2 + (L.b * cos(theta_ik[L.iO] - theta_ik[L.iD]) - L.g * sin(theta_ik[L.iO] - theta_ik[L.iD])) * v_ik[L.iO] * v_ik[L.iD] ))
            # 66
            push!(pD_e_constr_kl, @NLconstraint(OPF,
            pD_e_kl[end] == L.g * v_ik[L.iD]^2 + (-L.g * cos(theta_ik[L.iD] - theta_ik[L.iO]) - L.b * sin(theta_ik[L.iD] - theta_ik[L.iO])) * v_ik[L.iO] * v_ik[L.iD] ))
            # 67
            push!(qD_e_constr_kl, @NLconstraint(OPF,
            qD_e_kl[end] == -(L.b + L.bCH/2) * v_ik[L.iD]^2 + (L.b * cos(theta_ik[L.iD] - theta_ik[L.iO]) - L.g * sin(theta_ik[L.iD] - theta_ik[L.iO])) * v_ik[L.iO] * v_ik[L.iD] ))
        end

        sigS_e_kl = Array{JuMP.VariableRef, 1}()
        CurrentO_e_constr_kl = Array{JuMP.ConstraintRef, 1}()
        CurrentD_e_constr_kl = Array{JuMP.ConstraintRef, 1}()

        for L in LineList_k
            # 79
            push!(sigS_e_kl, @variable(OPF, lower_bound=0, base_name="sigS_$(L.e)_kl") )
            # 78
            push!(CurrentO_e_constr_kl, @NLconstraint(OPF, sqrt(pO_e_kl[end]^2 + qO_e_kl[end]^2) ≤ L.R_max*v_ik[L.iO] + sigS_e_kl[end]))
            # 80
            push!(CurrentD_e_constr_kl, @NLconstraint(OPF, sqrt(pD_e_kl[end]^2 + qD_e_kl[end]^2) ≤ L.R_max*v_ik[L.iD] + sigS_e_kl[end]))
        end
    end
    #
    # # there were no contigencies for transformers (but might be included in the future)
    #
    # # ∀ f ∈ F
    # # Transformer power ratings: 81-83
    # sigS_f_k = Array{JuMP.VariableRef, 1}()
    # PowerO_f_constr_k = Array{JuMP.ConstraintRef, 1}()
    # PowerD_f_constr_k = Array{JuMP.ConstraintRef, 1}()
    # for T in PN.TransformerList
    #     # 81(56)
    #     push!(sigS_f_k, @variable(OPF, lower_bound=0, base_name="sigS_$(T.f)_k") )
    #     # 82(55)
    #     push!(PowerO_f_constr_k, @NLconstraint(OPF, sqrt(pO_f_k[end]^2 + qO_f_k[end]^2) <= T.s_max + sigS_f_k[end]))
    #     # 83(57)
    #     push!(PowerD_f_const_k, @NLconstraint(OPF, sqrt(pD_f_k[end]^2 + qD_f_k[end]^2) <= T.s_max + sigS_f_k[end]))
    # end

    # # delta - powerfall in the case of contingency
    @variable(OPF, delta_k)
    @constraint(OPF, delta_k == sum(p_g_kg[i] for i = 1:sizeG) - sum(p_g[i] for i = 1:sizeG) )


    # Generator Real Power Contingency Response
    #(85)
    for i in 1:sizeG
        if i in PN.G
            push!(Pdelt_upper_constr, @constraint(OPF, (p_g_kg[end] - PN.GeneratorList[i].p_min)*(p_g_kg[end] - (p_g[end] + PN.GeneratorList[i].alpha*delta_k)) <= 0))
            push!(Pdelt_lower_constr, @constraint(OPF, (p_g_kg[end] - PN.GeneratorList[i].p_max)*(p_g_kg[end] - (p_g[end] + PN.GeneratorList[i].alpha*delta_k)) <= 0))
        end
    end

    # Generator Reactive Power Contingency Response
    #(93)
    for i in 1:sizeG
        if i in PN.G
            push!(QV_upper_constr, @constraint(OPF, (q_g_kg[end] - PN.GeneratorList[i].q_min)*(v_ik[end] - v_i[end]) <= 0))
            push!(QV_lower_constr, @constraint(OPF, (q_g_kg[end] - PN.GeneratorList[i].q_max)*(v_ik[end] - v_i[end]) <= 0))
        end
    end

    # optimize!(OPF)
    return OPF
end
