# model of the ARPA-E Competition
#
# JH: seems to be a bug somewhere.
#
# getting constraint violations for real and reactive power balance at all buses
# in the base case and all contingency cases.
#
# one possible problem is in the cos() and sin() expressions.
# I believe these expressions do not match the documentation exactly
# the documentation has the argument as
#
#   theta_to - theta_from + theta_transformer
#
# together with the +/- signs on the cos and sin terms that this model uses,
# I believe this model would be correct if the term
#
#   theta_transformer
#
# carried a factor of -1, i.e. a (-) sign rather than a (+) sign.
# not sure though as this depends on the data input process.
# Not correcting this yet.
# Better to see what the evaluation code comes up with after fixing
# the following error.
# it is possible that the code is correct, depending on the data input process.
#
# another possible problem is in the output of the voltage angles.
# it appears that the maximum angle difference in this model comes out as
# about 3 degrees (in solution2.txt), while the maximum angle difference in the
# gams model is 10 or 15 degrees.
# SOLVED: the issue here is the output to solution2.txt uses v and v0,
# rather than theta and theta0. This has been corrected.
#
# these observations are all without any pv/pq switching in the model.

using Complementarity;
using Ipopt;
using Gurobi;
#using KNITRO;
#using Ipopt, Complementarity;

# Global: contDList

function buildMod(fData,uData, contDList)
  # Input:
  # fData - grid data with all parameters
  # uData - contingency data

  # readin the data from the structure
  baseMVA = fData.baseMVA;
  bList = fData.busList;
  bData = fData.busDList;
#  println("bList:",bList); 	#STE

#  b0List = for i in bList if (bData[i].gsh !== 0)end end	#STE pFlow0[k in brList;brData[k].zeroImpe == false]
#  b0List = (bList[i] for i in bList if !(bData[i].gsh == 0))	#STE
#  b0List = for i in bList if bData[i].gsh == 0 end end	#STE
#  b0List = (i for i in bList if !(bData[i].gsh == 0))	#STE

#  b1List = (bList[i] for i in bList if (bData[i].gsh !== 0))	#STE

#  	println("b0List:",b0List); 	#STE
#  	println("b1List:",b1List); 	#STE

  gList = fData.genList;
  gData = fData.genDList;
 # 	println("gList:",gList); 	#STE

  brList = fData.brList;
  brData = fData.brDList;
#  	println("brList:",brList); 	#STE

  S = uData.contList;
  contData = uData.contDList;
#	println("S contList:",S); 	#STE

  # set up the model
  #mp = Model(solver = IpoptSolver(print_level=3));
  mp = Model(solver = IpoptSolver(max_iter=5000));
#mp = Model(solver = Gurobi());
  #mp = Model(solver = KnitroSolver(KTR_PARAM_OUTLEV=2,  # default is 2
 # 			bar_feasible=1,
 # 			bar_initmu=0.12,
 #  			feastol=2.25e-9,
   			#feastol_abs=1e-2,
 #  			opttol=1e-4,
 #  			cg_maxit=10,   # formerly maxcgit
 #  			maxit=400,
   			#ftol=1e-4,
   			#ftol_iters=3,
   			#pivot=1e-12,
 #  			maxtime_real=3600));

  # create the variables for the base case
 # bDV0 = @variable(mp,bData[i].Vmin <= v0[i in bList] <= bData[i].Vmax);	#STE 	Eqn. 19?
#  for i in bList; println("bData[",i,"].Vmin : ",bData[i].Vmin,", .Vmax: ",bData[i].Vmax); end;	#STE
#  println("bData[i].Vmin : ",bData[i].Vmin for i in bList);	#STE
#  println("bData[i].Vmin0: ",bDV0);	#STE

  gDP0 = @variable(mp,gData[l].Pmin <= sp0[l in gList] <= gData[l].Pmax);	#STE 	Eqn. 20?
#  for l in gList; println("gData[",l,"].Pmin : ",gData[l].Pmin,", .Pmax: ",gData[l].Pmax); end;	#STE
#  println("gData[l].Pmin0: ",gDP0);	#STE

  gDQ0 = @variable(mp,gData[l].Qmin <= sq0[l in gList] <= gData[l].Qmax);	#STE	Eqn. 21?
#  for l in gList; println("gData[",l,"].Qmin : ",gData[l].Qmin,", .Qmax: ",gData[l].Qmax); end;	#STE
#  println("gData[l].Qmin0: ",gDQ0);	#STE

  @variable(mp,p0[k in brList]);
#   println("mp.p0[k in brList]: ",mp);	#STE
  @variable(mp,q0[k in brList]);
#   println("mp.q0[k in brList]: ",mp);	#STE
  @variable(mp,psh0[i in bList]);
#   println("mp.psh0[i in bList]: ",mp);	#STE
  @variable(mp,qsh0[i in bList]);
#   println("mp.qsh0[i in bList]: ",mp);	#STE
#  @variable(mp,theta0[i in bList]);
#   println("mp.theta0[i in bList]: ",mp);	#STE
  @variable(mp,c0[i in bList, j in bList]);
  @variable(mp,s0[i in bList, j in bList]);

  # create the constraints for the base case
  # constraints for auxilary variables
  #@constraint(mp, symmetryc0[i in bList, j in bList], c0[i,j] == c0[j,i]);
  #@constraint(mp, symmetrys0[i in bList, j in bList], s0[i,j] == -s0[j,i]);
  # SOCP relaxation
  #@constraint(mp, SOCP[i in bList, j in bList], c0[i,j]^2 + s0[i,j]^2 + (0.5*(c0[i,i]-c0[j,j]))^2 ≤ (0.5*(c0[i,i]+c0[j,j]))^2);
  # Max 0.9 < = V <= 1.1
 # @constraint(mp, maxc0ii[i in bList], bData[i].Vmin^2 <= c0[i,i] <= bData[i].Vmax^2);

# @NLconstraint(mp,flowBound0[k in brList],p0[k]^2 + q0[k]^2 <=                 brData[k].t^2);				#STE original nonlinear inequality (1)
  @NLconstraint(mp,flowBound0[k in brList],p0[k]^2 + q0[k]^2 <= fData.baseMVA^2*brData[k].t^2);				#STE fixed    nonlinear inequality (1) Eqn. 43,44
# @NLconstraint(mp, pShunt0[   i in bList                      ], psh0[i] == bData[i].gsh * v0[i]^2);			#STE nonlinear original
#  @NLconstraint(mp, pShunt0_nl[i in bList; !(bData[i].gsh == 0)], psh0[i] == bData[i].gsh * v0[i]^2);			#STE nonlinear piece
  @NLconstraint(mp, pShunt0_nl[i in bList; !(bData[i].gsh == 0)], psh0[i] == bData[i].gsh * c0[i,i]);			#STE nonlinear piece
    @constraint(mp, pShunt0_li[i in bList;  (bData[i].gsh == 0)], psh0[i] == 0);					#STE    linear piece
# @NLconstraint(mp, qShunt0[   i in bList                      ], qsh0[i] == -bData[i].bsh * v0[i]^2);			#STE nonlinear original
#  @NLconstraint(mp, qShunt0_nl[i in bList; !(bData[i].bsh == 0)], qsh0[i] == -bData[i].bsh * v0[i]^2);			#STE nonlinear piece
@NLconstraint(mp, qShunt0_nl[i in bList; !(bData[i].bsh == 0)], qsh0[i] == -bData[i].bsh * c0[i,i]);			#STE nonlinear piece
    @constraint(mp, qShunt0_li[i in bList;  (bData[i].bsh == 0)], qsh0[i] == 0);  					#STE     linear piece
#=
  @NLconstraint(mp,pFlow0[k in brList;brData[k].zeroImpe == false], p0[k] == brData[k].g/(brData[k].tauprime^2)*v0[brData[k].From]^2
               - brData[k].g/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*cos(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               - brData[k].b/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*sin(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               );								#STE nonlinear
=#
@NLconstraint(mp,pFlow0[k in brList;brData[k].zeroImpe == false], p0[k] == brData[k].g/(brData[k].tauprime^2)*c0[brData[k].From,brData[k].From]
             - brData[k].g/brData[k].tau*c0[brData[k].From,brData[k].To]
             - brData[k].b/brData[k].tau*s0[brData[k].From,brData[k].To]
             );								#STE nonlinear
#=
  @NLconstraint(mp,qFlow0[k in brList;brData[k].zeroImpe == false], q0[k] == (-brData[k].b - brData[k].bc/2)/(brData[k].tauprime^2)*v0[brData[k].From]^2
               + brData[k].b/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*cos(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               - brData[k].g/brData[k].tau*v0[brData[k].From]*v0[brData[k].To]*sin(theta0[brData[k].From] - theta0[brData[k].To] + brData[k].thetatr)
               );								#STE nonlinear
=#
@NLconstraint(mp,qFlow0[k in brList;brData[k].zeroImpe == false], q0[k] == (-brData[k].b - brData[k].bc/2)/(brData[k].tauprime^2)*c0[brData[k].From,brData[k].From]
             + brData[k].b/brData[k].tau*c0[brData[k].From,brData[k].To]
             - brData[k].g/brData[k].tau*s0[brData[k].From,brData[k].To]
             );								#STE nonlinear
  @constraint(mp,zeroP0[k in brList;brData[k].zeroImpe == true], p0[k] == -p0[brData[k].revID]);	#STE linear equality
 # @NLconstraint(mp,zeroQ0[k in brList;brData[k].zeroImpe == true], q0[k] + q0[brData[k].revID] == -brData[k].bc*v0[brData[k].To]^2);
 @NLconstraint(mp,zeroQ0[k in brList;brData[k].zeroImpe == true], q0[k] + q0[brData[k].revID] == -brData[k].bc*c0[brData[k].To,brData[k].To]);
  #@constraint(mp,zeroV0[k in brList;brData[k].zeroImpe == true], v0[brData[k].To] == v0[brData[k].From]/brData[k].tau);
  #@constraint(mp,zerotheta0[k in brList;brData[k].zeroImpe == true], theta0[brData[k].To] == theta0[brData[k].From] - brData[k].thetatr);
  @constraint(mp,pBalance0[i in bList], sum(sp0[l] for l in bData[i].gen) == psh0[i] + bData[i].Pd + sum(p0[k] for k in brList if (brData[k].From == i)));
  @constraint(mp,qBalance0[i in bList], sum(sq0[l] for l in bData[i].gen) == qsh0[i] + bData[i].Qd + sum(q0[k] for k in brList if (brData[k].From == i)));
#    println("base constraints ",mp);	#STE

  # create the variables for the contingency cases
#  bDV = @variable(mp,bData[i].Vmin <= v[i in bList, s in S] <= bData[i].Vmax);
#  for i in bList; println("bData[",i,"].Vmin : ",bData[i].Vmin,", .Vmax: ",bData[i].Vmax); end;	#STE
#  println("bData[i].Vmin: ",bDV);	#STE

  gDP = @variable(mp,gData[l].Pmin <= sp[l in gList, s in S;(!(l in contData[s].Loc))] <= gData[l].Pmax);
#  for l in gList; println("gData[",l,"].Pmin : ",gData[l].Pmin,", .Pmax: ",gData[l].Pmax); end;	#STE
#  println("gData[l].Pmin: ",gDP);	#STE

  gDQ = @variable(mp,gData[l].Qmin <= sq[l in gList, s in S;(!(l in contData[s].Loc))] <= gData[l].Qmax);
#  for l in gList; println("gData[",l,"].Qmin : ",gData[l].Qmin,", .Qmax: ",gData[l].Qmax); end;	#STE
#  println("gData[l].Qmin: ",gDQ);	#STE

  @variable(mp,p[k in brList, s in S; !(k in contDList[s].Loc)]);
  @variable(mp,q[k in brList, s in S; !(k in contDList[s].Loc)]);
  @variable(mp,psh[i in bList, s in S]);
  @variable(mp,qsh[i in bList, s in S]);
 # @variable(mp,theta[i in bList, s in S]);
  @variable(mp,pdelta[s in S]);
  @variable(mp,ck[i in bList, j in bList, s in S]);
  @variable(mp,sk[i in bList, j in bList, s in S]);

  # create the constraints for the contingency cases
  # constraints for auxilary variables
 @constraint(mp, symmetryck[i in bList, j in bList, s in S], ck[i,j,s] == ck[j,i,s]);
 @constraint(mp, symmetrysk[i in bList, j in bList, s in S], sk[i,j,s] == -sk[j,i,s]);
  # SOCP relaxation
  @constraint(mp, SOCPk[i in bList, j in bList, s in S], ck[i,j,s]^2 + sk[i,j,s]^2 + (0.5*(ck[i,i,s]-ck[j,j,s]))^2 ≤ (0.5*(ck[i,i,s]+ck[j,j,s]))^2);
  #@constraint(mp, SOCPk[i in bList, j in bList, s in S], ck[i,j,s]^2 + sk[i,j,s]^2 == ck[i,i,s]*ck[j,j,s]);
  # Max 0.9 < = V <= 1.1
  @constraint(mp, maxckii[i in bList, s in S], bData[i].Vmin^2 <= ck[i,i,s] <= bData[i].Vmax^2);
  # Phase angle constraint: θi - θj <=10°
  @constraint(mp, LimitCij[i in bList, j in bList, s in S], 0.79 ≤ ck[i,j,s] ≤ 1.21);
  @constraint(mp, LimitSij[i in bList, j in bList, s in S], -0.21 ≤ sk[i,j,s] ≤ 0.21);
  # 3 cycles
b_keys_sorted = sort(collect(keys(bData)));
NBus = size(bList,1);
NBr = size(brList,1);
# Bus-bus incidence matrix
Ψ = zeros(NBus,NBus);
for i = 1:NBus
    for j = 1:NBr
        if collect(bData[k].ID for k in b_keys_sorted)[i] == brList[j][1]
            Ψ[i,brList[j][2]] = 1
        end
    end
end
Ψh = UpperTriangular(Ψ);

cycle3list = [];
for i in 1:NBus
    for j in 1:NBus
        for k in 1:NBus
            if Ψh[i,j] == 1 && Ψh[i,k] == 1 && Ψh[k,j] == 1
                cycle3 = [i,j,k]
                push!(cycle3list,cycle3)
            end
        end
    end
end
cycle3_size  = size(cycle3list)[1];

#Original q13 and q23 constraints
@constraint(mp, q13[i = 1:cycle3_size, s in S], sk[cycle3list[i][1],cycle3list[i][2],s] * ck[cycle3list[i][3],cycle3list[i][3],s]
+ ck[cycle3list[i][2],cycle3list[i][3],s] * sk[cycle3list[i][3],cycle3list[i][1],s]
+ sk[cycle3list[i][2],cycle3list[i][3],s] * ck[cycle3list[i][3],cycle3list[i][1],s] == 0);

@constraint(mp, q23[i = 1:cycle3_size, s in S], ck[cycle3list[i][1],cycle3list[i][2],s] * ck[cycle3list[i][3],cycle3list[i][3],s]
- ck[cycle3list[i][2],cycle3list[i][3],s] * ck[cycle3list[i][3],cycle3list[i][1],s]
+ sk[cycle3list[i][2],cycle3list[i][3],s] * sk[cycle3list[i][3],cycle3list[i][1],s] == 0);

# @NLconstraint(mp,flowBoundS[k in brList, s in S; !(k in contDList[s].Loc)],p[k,s]^2 + q[k,s]^2 <=                 brData[k].t^2);	#STE original
  @NLconstraint(mp,flowBoundS[k in brList, s in S; !(k in contDList[s].Loc)],p[k,s]^2 + q[k,s]^2 <= fData.baseMVA^2*brData[k].t^2);	#STE fixed
# @NLconstraint(mp, pShuntS[   i in bList, s in S]                      , psh[i,s] == bData[i].gsh * v[i,s]^2); 	#STE original
#  @NLconstraint(mp, pShuntS_nl[i in bList, s in S; !(bData[i].gsh == 0)], psh[i,s] == bData[i].gsh * v[i,s]^2);		#STE nonlinear piece
  @NLconstraint(mp, pShuntS_nl[i in bList, s in S; !(bData[i].gsh == 0)], psh[i,s] == bData[i].gsh * ck[i,i,s]);
    @constraint(mp, pShuntS_li[i in bList, s in S;  (bData[i].gsh == 0)], psh[i,s] == 0);				#STE    linear piece
# @NLconstraint(mp, qShuntS[   i in bList, s in S],			  qsh[i,s] == -bData[i].bsh * v[i,s]^2);	#STE original
  #@NLconstraint(mp, qShuntS_nl[i in bList, s in S; !(bData[i].bsh == 0)], qsh[i,s] == -bData[i].bsh * v[i,s]^2);	#STE nonlinear piece
  @NLconstraint(mp, qShuntS_nl[i in bList, s in S; !(bData[i].bsh == 0)], qsh[i,s] == -bData[i].bsh * ck[i,i,s]);	#STE nonlinear piece
    @constraint(mp, qShuntS_li[i in bList, s in S;  (bData[i].bsh == 0)], qsh[i,s] == 0);				#STE    linear piece
#=
  @NLconstraint(mp,pFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], p[k,s] == brData[k].g/(brData[k].tauprime^2)*v[brData[k].From,s]^2
               - brData[k].g/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*cos(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               - brData[k].b/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*sin(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               );
=#
    @NLconstraint(mp,pFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], p[k,s] == brData[k].g/(brData[k].tauprime^2)*ck[brData[k].From,brData[k].From,s]
                - brData[k].g/brData[k].tau*ck[brData[k].From,brData[k].To,s]
                - brData[k].b/brData[k].tau*sk[brData[k].From,brData[k].To,s]
                );
#=
  @NLconstraint(mp,qFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], q[k,s] == (-brData[k].b - brData[k].bc/2)/(brData[k].tauprime^2)*v[brData[k].From,s]^2
               + brData[k].b/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*cos(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               - brData[k].g/brData[k].tau*v[brData[k].From,s]*v[brData[k].To,s]*sin(theta[brData[k].From,s] - theta[brData[k].To,s] + brData[k].thetatr)
               );
=#
@NLconstraint(mp,qFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], q[k,s] == (-brData[k].b - brData[k].bc/2)/(brData[k].tauprime^2)*ck[brData[k].From,brData[k].From,s]
             + brData[k].b/brData[k].tau*ck[brData[k].From,brData[k].To,s]
             - brData[k].g/brData[k].tau*sk[brData[k].From,brData[k].To,s]
             );
  zeroPSC = @constraint(mp,zeroPS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], p[k,s] == -p[brData[k].revID,s]);
  #zeroQSC = @NLconstraint(mp,zeroQS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], q[k,s] + q[brData[k].revID,s] == -brData[k].bc*v[brData[k].To,s]^2);
  zeroQSC = @NLconstraint(mp,zeroQS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], q[k,s] + q[brData[k].revID,s] == -brData[k].bc*ck[brData[k].To,brData[k].To,s]);
  #@constraint(mp,zeroVS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], v[brData[k].To,s] == v[brData[k].From,s]/brData[k].tau);
  #@constraint(mp,zerothetaS[k in brList, s in S;(brData[k].zeroImpe == true)&(!(k in contDList[s].Loc))], theta[brData[k].To,s] == theta[brData[k].From,s] - brData[k].thetatr);
  @constraint(mp,pBalanceS[i in bList, s in S], sum(sp[l,s] for l in bData[i].gen if (!(l in contData[s].Loc))) == psh[i,s] + bData[i].Pd + sum(p[k,s] for k in brList if ((brData[k].From == i)&(!(k in contDList[s].Loc)))));
  @constraint(mp,qBalanceS[i in bList, s in S], sum(sq[l,s] for l in bData[i].gen if (!(l in contData[s].Loc))) == qsh[i,s] + bData[i].Qd + sum(q[k,s] for k in brList if ((brData[k].From == i)&(!(k in contDList[s].Loc)))));
  recourseC = @constraint(mp,recourse[l in gList, s in S; (!(l in contData[s].Loc))], sp[l,s] == sp0[l] + (gData[l].alpha)*pdelta[s]);	#STE! linear equality
#  println("contingency constraints: zeroPS ",zeroPSC); #STE
#  println("contingency constraints: zeroQS ",zeroQSC); #STE
#  println("contingency constraints: recourse() ",recourseC); #STE
#    println("contingency constraints ",mp);	#STE

  # create the Complementarity constraint
  @variable(mp, delplus[i in bList, s in S] >= 0);
  @variable(mp, delminus[i in bList, s in S] >= 0);
  @variable(mp, sqplus[i in bList,s in S] >= 0);
  @variable(mp, sqminus[i in bList,s in S] >= 0);

  #@constraint(mp,vConstr1[i in bList,s in S], v0[i] - v[i,s] <= delplus[i,s]);
  #@constraint(mp,vConstr2[i in bList,s in S], v[i,s] - v0[i] <= delminus[i,s]);
  @constraint(mp, spplusConstr[i in bList,s in S], sqplus[i,s] == sum(gData[l].Qmax - sq[l,s] for l in bData[i].gen));
  @constraint(mp, spminusConstr[i in bList,s in S], sqminus[i,s] == sum(sq[l,s] - gData[l].Qmin for l in bData[i].gen));
#    println("complementary constraints ",mp);	#STE
#=
  for s in S
    for i in bList
      @complements(mp,0 <= delplus[i,s], sqplus[i,s] >= 0);
      @complements(mp,0 <= delminus[i,s], sqminus[i,s] >= 0);
#      println("s ",s," i ",i," delplus[i,s] ",delplus[i,s]," delminus[i,s] ",delminus[i,s]," sqplus[i,s] ",sqplus[i,s]," sqminus[i,s] ",sqminus[i,s])
    end
  end
=#
#    println("complements ",mp);	#STE
#  println("complementarity is off");
#  # add the contingency effect as constraints
#  for s in S
#    # for each contingency
#    # first identify the type of the contingency
#    if contData[s].Type == "B"
#      # if it is a branch contingency
#      for k in contData[s].Loc
#        @constraint(mp, p[k,s] == 0.0);       # the real power flow through the branch is 0
#        @constraint(mp, q[k,s] == 0.0);       # the reactive power flow through the branch is 0
#      end
#    elseif contData[s].Type == "T"
#      # if it is a branch contingency
#      for k in contData[s].Loc
#        @constraint(mp, p[k,s] == 0.0);       # the real power flow through the branch is 0
#        @constraint(mp, q[k,s] == 0.0);       # the reactive power flow through the branch is 0
#      end
#    elseif contData[s].Type == "G"
#      # if it is a generator contingency
#      for l in contData[s].Loc
#        @constraint(mp, sp[l,s] == 0.0);
#        @constraint(mp, sq[l,s] == 0.0);
#      end
#    end
#  end

# 	build the objective function
#
# @variable(mp, CF); #STE
#  @NLconstraint(mp, CF == sum(sum(gData[l].cParams[n]*(sp0[l]*fData.baseMVA)^n for n in gData[l].cn) for l in gList if gData[l].Pmin<gData[l].Pmax)			#STE
#                        + sum(sum(gData[l].cParams[n]*(gData[l].Pmin*fData.baseMVA)^n for n in gData[l].cn ) for l in gList if gData[l].Pmin==gData[l].Pmax) );	#STE
#  @NLobjective(mp, Min, CF);  #STE needed to avoid error for RTS96
@NLobjective(mp, Min,    sum(sum(gData[l].cParams[n]*(sp0[l]*fData.baseMVA)^n for n in gData[l].cn) for l in gList if gData[l].Pmin<gData[l].Pmax)			#STE
                         + sum(sum(gData[l].cParams[n]*(gData[l].Pmin*fData.baseMVA)^n for n in gData[l].cn ) for l in gList if gData[l].Pmin==gData[l].Pmax) );	#STE
# @NLobjective(mp, Min,    sum(sum(gData[l].cParams[n]*(sp0[l]*fData.baseMVA)^n for n in gData[l].cn) for l in gList)); #original

  return mp;
end
