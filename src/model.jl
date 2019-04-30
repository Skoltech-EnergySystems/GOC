using JuMP
using Ipopt
include("NetworkData.jl")


function create_model(PN::PNetwork)

    # OPF = Model(solver=IpoptSolver())
    OPF = Model(with_optimizer(Ipopt.Optimizer))


    # base model

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

    for g in PN.G
        i = PN.gen_ind_I[g] # index of generator[g] in the list of all generators
    # when i wtire this constraint i set some limits like 2 <= p_g <= 10
    # and then i add constraints for inactive generators like p_g == 0 which is a constradiction
    # with the constrain above.
# 33
        push!(p_g, @variable(OPF, lower_bound=PN.GeneratorList[i].p_min, upper_bound=PN.GeneratorList[i].p_max, base_name="p_$g"))
# 35
        push!(q_g, @variable(OPF, lower_bound=PN.GeneratorList[i].q_min, upper_bound=PN.GeneratorList[i].q_max, base_name="q_$g"))
# 4, 5
        push!(t_gh, @variable(OPF, [h=1:PN.GeneratorList[i].N], lower_bound=0, base_name="t_$g"))
        push!(t_g_constr, @constraint(OPF, sum(t_gh[end]) == 1))
# 34, 36
        if g in diffG
            push!(p_g_constr, @constraint(OPF, p_g[end] == 0))
            push!(q_g_constr, @constraint(OPF, q_g[end] == 0))
        else
# 3
            push!(p_g_constr, @constraint(OPF, p_g[end] == PN.GeneratorList[i].p' * t_gh[end]))
        end
# 2
        push!(c_g, @variable(OPF, base_name="c_$g"))
        push!(c_g_constr, @constraint(OPF, c_g[end] == PN.GeneratorList[i].c' * t_gh[end]))
    end

# ∀ i ∈ I
    v_i = Array{JuMP.VariableRef, 1}()
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
        if i in keys(PN.get_shunt_index)
            push!(bCS_i, @variable(OPF, lower_bound=PN.sShuntList[PN.get_shunt_index[i]].bCS_min, upper_bound=PN.sShuntList[PN.get_shunt_index[i]].bCS_max, base_name="bCS_$(PN.get_shunt_index[i])"))
        end
# 47 48 50 51
        push!(sigPp_i, @variable(OPF, lower_bound=0, base_name="sigPp_$i") )
        push!(sigPm_i, @variable(OPF, lower_bound=0, base_name="sigPm_$i") )
        push!(sigQp_i, @variable(OPF, lower_bound=0, base_name="sigQp_$i") )
        push!(sigQm_i, @variable(OPF, lower_bound=0, base_name="sigQm_$i") )
# 46
        # push!(sigP_diff, @constraint(OPF, sigPp_i[end] - sigPm_i[end] == )  )

    end


size = maximum(PN.caliI);
@variable(OPF, theta[1:size])
@variable(OPF, v[1:size])
pO_e = Array{JuMP.VariableRef, 1}()
pO_constr = Array{JuMP.ConstraintRef, 1}()
# ∀ e ∈ E
for L in PN.LineList
    push!(pO_e, @variable(OPF, base_name="pO_e_$(L.e)"))
# 38
    push!(pO_constr, @NLconstraint(OPF,
    pO_e[end] == L.g * v[L.iO]^2 + (-L.g * cos(theta[L.iO] - theta[L.iD]) - L.b * sin(theta[L.iO] - theta[L.iD])) * v[L.iO] * v[L.iD] ))
end


# ∀ f ∈ F





    # model for contingencies

    return OPF
end
