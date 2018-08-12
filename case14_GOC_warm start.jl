using PowerModels
using Ipopt
using CSV
solver = IpoptSolver()
cd("C:/Users/Ильгиз/Documents/Документы_Ильгиз/Skoltech_2018/Grid Competition/Code")
############################ GOC DATA k1 #####################################
result1 = run_ac_opf("case14_GOC_6-12 off.m",solver)
# Objective function value
#println("RESULTS: ", result1)
println("objective function value: ",result1["objective"])

# Optimal voltage magnitude and phase angle
bus_info1 = sort(result1["solution"]["bus"]);
println("Optimal voltage phase angle and magnitude")
for (key, value) in bus_info1
    println(key, " ==> ", value)
 end

# Optimal generation
gen_info1 = sort(result1["solution"]["gen"]);
println("Optimal reactive and active generation")
for (key, value) in gen_info1
    println(key, " ==> ", value)
 end
############################ GOC DATA k0 #####################################
result2 = run_ac_opf("case14_GOC.m",solver)
# Objective function value
println("objective function value: ",result2["objective"])

# Optimal voltage magnitude and phase angle
bus_info2 = sort(result2["solution"]["bus"]);
println("Optimal voltage phase angle and magnitude")
for (key, value) in bus_info2
    println(key, " ==> ", value)
 end

# Optimal generation
gen_info2 = sort(result2["solution"]["gen"]);
println("Optimal reactive and active generation")
for (key, value) in gen_info2
    println(key, " ==> ", value)
 end
