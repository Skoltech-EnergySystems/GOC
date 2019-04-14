# Writing into the file "solution1.txt"

NBus = size(Bus_I,1)
sSh_ind = [s.i for s in PN.sShuntList]
NGen = size(Gen_I,1)

Nod = zeros(NBus,4)
Nod[1:NBus,1] = Bus_I[:,1]
Nod[1:NBus,2] = V_I[:,1]
Nod[1:NBus,3] = Theta_I[:,1]
Nod[sSh_ind,4] = BCS_I

Gen = zeros(NGen,4)
Gen[1:NGen,1] =Gen_I[:,1]
Gen[1:NGen,2] =Gen_ID[:,1]
Gen[1:NGen,3] =P_G[:,1]
Gen[1:NGen,4] =Q_G[:,1]

open("solution1.txt", "w") do f1
    write(f1, "--bus section\r\n")
    write(f1, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\r\n")
    for i in 1:NBus

        write(f1, join(Nod[i,:], ","))
        write(f1, "\r\n")
    end
    write(f1, "--generator section\r\n")
    write(f1, "i, id, p(MW), q(MVAR)\r\n")
    for i in 1:NGen
        write(f1, join(Gen[i,:], ","))
        write(f1, "\r\n")
    end
end
