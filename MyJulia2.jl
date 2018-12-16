include("def.jl");
include("dProc.jl");
include("buildMod.jl");

function MyJulia2(rawFile, genFile, contFile) 

# written by Haoxiang Yang, Northwestern University, 404-421-0638, haoxiangyang2019@u.northwestern.edu
# modified by Stephen Elbert and Jesse Holzer, PNNL; Miles Lubin, Google
   
    println("RAW file: ",rawFile)
    println("gen file: ",genFile)
    println("con file: ",contFile)

    println("Data generation")

    tic()

    baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg = preProc(rawFile);
    busList,busDList,genList,genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
    brList,brListSingle,brDList = branchProc(baseMVA,branchSeg,transformerSeg);
    contList,contDList = contProc(contFile);

    fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
    uData = uncertainData(contList,contDList);

    toc()

    println("Model generation")

    tic()
    mp = buildMod(fData,uData, contDList);
    toc()

    println("Model solve")
    tic()
    solve(mp);
    toc()
    sphat = getvalue(mp[:sp0]);
    sqhat = getvalue(mp[:sq0]);
    spshat = getvalue(mp[:sp]);
    sqshat = getvalue(mp[:sq]);
    costhat = getobjectivevalue(mp);

    println("Solution writing")
    tic()

    open("solution1.txt","w") do f
      write(f, "--generation dispatch \nbus id,unit id,pg(MW),qg(MVar) \n");
      for i in fData.genList
        loc = fData.genDList[i].Loc;
        name = fData.genDList[i].Name;
        spTemp = sphat[i]*fData.baseMVA;
        sqTemp = sqhat[i]*fData.baseMVA;
        write(f, "$loc,$name,$spTemp,$sqTemp \n");
      end
      write(f,"--end of generation dispatch \n");
    end

    open("solution2.txt","w") do f
      write(f, "--contingency generator \nconID,genID,busID,unitID,q(MW) \n");
      for s in uData.contList
        counter = 0;
        for i in fData.genList
          counter += 1;
          loc = fData.genDList[i].Loc;
          idTemp = fData.genDList[i].ID;
          name = fData.genDList[i].Name;
          sqTemp = sqshat[i,s]*fData.baseMVA;
          write(f,"$s,l_$counter,$loc,$name,$sqTemp \n");
        end
      end
      write(f,"--end of contingency generator \n--bus \ncontingency id,bus id,v(pu),theta(deg) \n");
      for i in fData.busList
        id = fData.busDList[i].ID;
        name = fData.busDList[i].Name;
        vTemp = getvalue(mp[:v0][i]);
        thetaTemp = getvalue(mp[:theta0][i]/pi*180);
        write(f,"0,$id,$vTemp,$thetaTemp \n");
      end
      for s in uData.contList
        for i in fData.busList
          id = fData.busDList[i].ID;
          name = fData.busDList[i].Name;
          vTemp = getvalue(mp[:v][i,s]);
          thetaTemp = getvalue(mp[:theta][i,s]/pi*180);
          write(f,"$s,$id,$vTemp,$thetaTemp \n");
        end
      end
      write(f,"--end of bus \n--Delta \ncontingency id,Delta(MW) \n");
      for s in uData.contList
        pdeltaTemp = getvalue(mp[:pdelta][s])*fData.baseMVA;
        write(f,"$s,$pdeltaTemp \n");
      end
      write(f,"--end of Delta \n--line flow \ncontingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar) \n");
      for i in fData.brListSingle
        idTemp = fData.brDList[i].ID;
        revidTemp = fData.brDList[i].revID;
        fromTemp = fData.brDList[i].From;
        toTemp = fData.brDList[i].To;
        name = fData.brDList[i].CKT;
        pTemp = getvalue(mp[:p0][i])*fData.baseMVA;
        qTemp = getvalue(mp[:q0][i])*fData.baseMVA;
        pRevTemp = getvalue(mp[:p0][revidTemp])*fData.baseMVA;
        qRevTemp = getvalue(mp[:q0][revidTemp])*fData.baseMVA;
        write(f,"0,$name,$fromTemp,$toTemp,i_$(fromTemp)_$(toTemp)_$(name),$pTemp,$qTemp,$pRevTemp,$qRevTemp \n");
      end
      for s in uData.contList
        for i in fData.brListSingle
          if !(i in contDList[s].Loc)
            idTemp = fData.brDList[i].ID;
            revidTemp = fData.brDList[i].revID;
            fromTemp = fData.brDList[i].From;
            toTemp = fData.brDList[i].To;
            name = fData.brDList[i].CKT;
            pTemp = getvalue(mp[:p][i,s])*fData.baseMVA;
            qTemp = getvalue(mp[:q][i,s])*fData.baseMVA;
            pRevTemp = getvalue(mp[:p][revidTemp,s])*fData.baseMVA;
            qRevTemp = getvalue(mp[:q][revidTemp,s])*fData.baseMVA;
            write(f,"$s,$name,$fromTemp,$toTemp,i_$(fromTemp)_$(toTemp)_$(name),$pTemp,$qTemp,$pRevTemp,$qRevTemp \n");
          end
        end
      end
      write(f,"--end of line flow \n")
    end

    toc()
end
