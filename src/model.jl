using JuMP
using Ipopt
include("NetworkData.jl")


function create_model(PN::PNetwork)

    OPF = Model(solver=IpoptSolver())

    # base model

    # """
    # loop through all complex indeces like g = (i, id) and add constraints manually.
    # """
    # base
    # @variables OPF begin
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
    # @variables OPF begin
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

    # variables and constraints for generators ∀ g ∈ G
    p_g = Array{JuMP.Variable, 1}()
    q_g = Array{JuMP.Variable, 1}()
    t_gh = Array{Array{JuMP.Variable, 1}, 1}()
    c_g = Array{JuMP.Variable, 1}()

    lenG = length(PN.G)
    @constraintref p_g_constr[1:lenG]
    @constraintref q_g_constr[1:lenG]
    @constraintref c_g_constr[1:lenG]
    @constraintref t_g_constr[1:lenG]

    diffG = setdiff(PN.caliG, PN.G)
    j = 1
    for g in PN.G
        i = PN.gen_ind_I[g] # index of generator[g] in the list of all generators

        push!(p_g, @variable(OPF, lowerbound=PN.GeneratorList[i].p_min, upperbound=PN.GeneratorList[i].p_max, basename="p_$g"))

        push!(q_g, @variable(OPF, lowerbound=PN.GeneratorList[i].q_min, upperbound=PN.GeneratorList[i].q_max, basename="q_$g"))

        push!(t_gh, @variable(OPF, [h=1:PN.GeneratorList[i].N], lowerbound=0, basename="t_$g"))
        # constraint on sum of t_gh for each g
        t_g_constr[j] = @constraint(OPF, sum(t_gh[end]) == 1)

        if g in diffG
            p_g_constr[j] = @constraint(OPF, p_g[end] == 0)
            q_g_constr[j] = @constraint(OPF, q_g[end] == 0)
        else
            p_g_constr[j] = @constraint(OPF, p_g[end] == PN.GeneratorList[i].p' * t_gh[end])
        end
        push!(c_g, @variable(OPF, basename="c_$g"))
        c_g_constr[j] = @constraint(OPF, c_g[end] == PN.GeneratorList[i].c' * t_gh[end])

        j += 1
    end


    # ∀ i ∈ I
    v_i = Array{JuMP.Variable, 1}()
    bCS_i = Array{JuMP.Variable, 1}()


    j = 1
    for i in PN.caliI
        push!(v_i, @variable(OPF, lowerbound=PN.BusList[i].v_min, upperbound=PN.BusList[i].v_max, basename="v_$i"))

        if i in keys(PN.get_shunt_index)
            push!(bCS_i, @variable(OPF, lowerbound=PN.sShuntList[PN.get_shunt_index[i]].bCS_min, upperbound=PN.sShuntList[PN.get_shunt_index[i]].bCS_max, basename="bCS_$(PN.get_shunt_index[i])"))
        end

        j += 1
    end






    # model for contingencies

    return OPF
end
