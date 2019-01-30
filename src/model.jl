using JuMP
using Ipopt
include("NetworkData.jl")

function init_model(PN::PNetwork, PC)

OPF = Model(solver=IpoptSolver())

# base model

# """
# loop through all complex indeces like g = (i, id) and add constraints manually.
# """
# base
@variables OPF begin
    c
    c_g
    cSigma
    t_gh
    v_i
    theta_i
    bCS_i
    sigPp_i
    sigPm_i
    sigQp_i
    sigQm_i
    sigS_e
    sigS_f
    p_g
    q_g
    pO_e
    pD_e
    qO_e
    qD_e
    pO_f
    pD_f
    qO_f
    qD_f
end

# contingency
@variables OPF begin
    cSig_k
    delta_k
    v_ik
    theta_ik
    bCS_ik
    sigPp_ik
    sigPm_ik
    sigQp_ik
    sigQm_ik
    sigS_ek
    sigS_fk
    p_gk
    q_gk
    pO_ek
    pD_ek
    qO_ek
    qD_ek
    pO_fk
    pD_fk
    qO_fk
    qD_fk
end


for g in PN.G
    
end




# model for contingencies
end
