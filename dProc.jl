# data processing for ARPA-E Competition Phase 0

# preprocessing: read in the data and divide it into segments
function preProc(rawFile)
  rawData = readdlm(rawFile,',', skipblanks=false);

  # separate the data into different segments
  titleLine1 = rawData[1,:];
#  println("titleLine1: ",titleLine1);
  titleLine2 = rawData[2,:];
#  println("titleLine2: ",titleLine2);
  titleLine3 = rawData[3,:];
#  println("titleLine3: ",titleLine3);
  n,m = size(rawData);
  busStartL = 4; # first three lines are headers
  busEndL = 0
  loadStartL = 0
  loadEndL = 0
  shuntStartL = 0
  shuntEndL = 0
  genStartL = 0
  genEndL = 0
  branchStartL = 0
  branchEndL = 0
  transformerStartL = 0
  transformerEndL = 0
  for i in 1:n
    if contains(string(rawData[i,1]),"END OF BUS DATA")
      busEndL = i - 1;
      loadStartL = i + 1;
    end
    if contains(string(rawData[i,1]),"END OF LOAD DATA")
      loadEndL = i - 1;
      shuntStartL = i + 1;
    end
    if contains(string(rawData[i,1]),"END OF FIXED SHUNT DATA")
      shuntEndL = i - 1;
      genStartL = i + 1;
    end
    if contains(string(rawData[i,1]),"END OF GENERATOR DATA")
      genEndL = i - 1;
      branchStartL = i + 1;
    end
    if contains(string(rawData[i,1]),"END OF BRANCH DATA")
      branchEndL = i - 1;
      transformerStartL = i + 1;
    end
    if contains(string(rawData[i,1]),"END OF TRANSFORMER DATA")
      transformerEndL = i - 1;
    end
  end

  baseMVA = rawData[1,2];
  busSeg = rawData[busStartL:busEndL,:];
#  println("busSeg: ",busSeg);
  loadSeg = rawData[loadStartL:loadEndL,:];
  shuntSeg = rawData[shuntStartL:shuntEndL,:];
  genSeg = rawData[genStartL:genEndL,:];
  branchSeg = rawData[branchStartL:branchEndL,:];
  transformerSeg = rawData[transformerStartL:transformerEndL,:];

  return baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg;
end

# parse the data for each bus and generator
function busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile)
  # process the busSeg
  bn,bm = size(busSeg);
  gn,gm = size(genSeg);
  sn,sm = size(shuntSeg);
  ln,lm = size(loadSeg);

  busList = [];
  busDList = Dict();
  busNo = 0;
  for i in 1:bn
    # build the bus item
    busID = busSeg[i,1];
#  println("busID: ",busID);
    push!(busList,busID);
    busName = string(busSeg[i,2]);
#  println("busName: ",busName);
    busVmax = busSeg[i,10];
#  println("busVmax: ",busVmax);
    busVmin = busSeg[i,11];
#  println("busVmin: ",busVmin);
    busItem = busData(busID,busName,[],busVmax,busVmin,0,0,0,0);
#  println("busItem: ",busItem);

    # add the bus item to the dictionary
    busDList[busID] = busItem;
  end

  # process the shuntSeg
  for i in 1:sn
    busID = shuntSeg[i,1];
    busDList[busID].gsh += shuntSeg[i,4]/baseMVA;
    busDList[busID].bsh += shuntSeg[i,5]/baseMVA;
  end

  # process the loadSeg
  for i in 1:ln
    busID = loadSeg[i,1];
    busDList[busID].Pd += loadSeg[i,6]/baseMVA;
    busDList[busID].Qd += loadSeg[i,7]/baseMVA;
  end

  # process the cost data from genFile
  genCostData,genCostTitle = readdlm(genFile,',',header = true);
  gcn,gcm = size(genCostData);
  genCost = Dict();
#  println("gcn: ",gcn);
#  println("gcm: ",gcm);
#  println("genCostData: ",genCostData);
  for i in 1:gcn
    if !((genCostData[i,1],genCostData[i,2]) in keys(genCost))
      genCost[(genCostData[i,1],genCostData[i,2])] = Dict();
      genCost[(genCostData[i,1],genCostData[i,2])][genCostData[i,3]] = genCostData[i,4];
    else
      genCost[(genCostData[i,1],genCostData[i,2])][genCostData[i,3]] = genCostData[i,4];
    end
  end

  # process the genSeg
  genList = [];
  genDList = Dict();
  for i in 1:gn
    busID = genSeg[i,1];
    genLoc = busID;
    genName = genSeg[i,2];
    genID = (genLoc,genName);
    push!(genList,genID);
    push!(busDList[busID].gen,genID);

    genQmax = genSeg[i,5]/baseMVA;
    genQmin = genSeg[i,6]/baseMVA;
    genPmax = genSeg[i,17]/baseMVA;
    genPmin = genSeg[i,18]/baseMVA;

    genCn = [];
    gencParams = Dict();
    genalpha = 1;
    for j in keys(genCost[(genLoc,genName)])
      if j != 9
        push!(genCn,j);
        gencParams[j] = genCost[(genLoc,genName)][j];
      else
        genalpha = genCost[(genLoc,genName)][j];
      end
    end

    genItem = genData(genID,genName,genLoc,genCn,gencParams,genPmax,genPmin,genQmax,genQmin,genalpha);
    genDList[genID] = genItem;
  end

  return busList,busDList,genList,genDList;
end
#println(busDList)
# parse the data for each branch
function branchProc(baseMVA,branchSeg,transformerSeg)
  # process the branchSeg
  brn,brm = size(branchSeg);
  brList = [];
  brListSingle = [];
  brDList = Dict();
  for i in 1:brn
    brFrom = branchSeg[i,1];
    brTo = branchSeg[i,2];
    brCKT = branchSeg[i,3];
    branchID = (brFrom,brTo,brCKT);
    bRevID = (brTo,brFrom,brCKT);

    brr = branchSeg[i,4];
    brx = branchSeg[i,5];
    brg = brr/(brr^2 + brx^2);
    brb = -brx/(brr^2 + brx^2);
    brbc = branchSeg[i,6];
    brt = branchSeg[i,7]/baseMVA;

    if (brr == 0)&(brx == 0)
      brZeroImpe = true;
    else
      brZeroImpe = false;
    end

    brtau = 1;
    brtauprime = 1;
    brthetatr = 0;
    brItem = branchData(brFrom,brTo,brCKT,branchID,bRevID,brr,brx,brg,brb,brbc,brt,brZeroImpe,brtau,brtauprime,brthetatr);
    brRevItem = branchData(brTo,brFrom,brCKT,bRevID,branchID,brr,brx,brg,brb,brbc,brt,brZeroImpe,brtau,brtauprime,brthetatr);
    push!(brList,branchID);
    push!(brListSingle,branchID);
    push!(brList,bRevID);
    brDList[branchID] = brItem;
    brDList[bRevID] = brRevItem;
  end

  # process the transformerSeg
  trn,trm = size(transformerSeg);
  lineNo = 1;
  while lineNo <= trn
    trFrom = transformerSeg[lineNo,1];
    trTo = transformerSeg[lineNo,2];
    if transformerSeg[lineNo,3] == 0
      lineNoNew = lineNo + 4;
    end
    trName = transformerSeg[lineNo,4];
    trID = (trFrom,trTo,trName);
    trRevID = (trTo,trFrom,trName);

    trr = transformerSeg[lineNo+1,1];
    trx = transformerSeg[lineNo+1,2];
    trg = trr/(trr^2 + trx^2);
    trb = -trx/(trr^2 + trx^2);
    trbc = 0;
    trt = transformerSeg[lineNo+2,4]/baseMVA;

    if (trr == 0)&(trx == 0)
      trZeroImpe = true;
    else
      trZeroImpe = false;
    end

    trthetatr = -transformerSeg[lineNo+2,3];
    trRevthetatr = transformerSeg[lineNo+2,3];
    trtau = transformerSeg[lineNo+2,1]/transformerSeg[lineNo+3,1];
    trtauprime = transformerSeg[lineNo+2,1]/transformerSeg[lineNo+3,1];
    trRevtau = transformerSeg[lineNo+2,1]/transformerSeg[lineNo+3,1];
    trRevtauprime = 1;
    trItem = branchData(trFrom,trTo,trName,trID,trRevID,trr,trx,trg,trb,trbc,trt,trZeroImpe,trtau,trtauprime,trthetatr);
    trRevItem = branchData(trTo,trFrom,trName,trRevID,trID,trr,trx,trg,trb,trbc,trt,trZeroImpe,trRevtau,trRevtauprime,trRevthetatr);

    push!(brList,trID);
    push!(brListSingle,trID);
    push!(brList,trRevID);
    brDList[trID] = trItem;
    brDList[trRevID] = trRevItem;
    lineNo = lineNoNew;
  end

  return brList,brListSingle,brDList;
end

# process the contingency data
function contProc(contFile)
  contData,contTitle = readdlm(contFile,',',header = true);
  contn,contm = size(contData);
  contList = [];
  contDList = Dict();
  for i in 1:contn
    contID = contData[i,1];
    push!(contList,contID);
    contType = contData[i,2];
    if (contType == "B")|(contType == "T")
      contLoc = [(contData[i,3],contData[i,4],contData[i,5]),(contData[i,4],contData[i,3],contData[i,5])];
    else
      contLoc = [contData[i,3],contData[i,4]];
    end
    contItem = contingencyData(contID,contType,contLoc);
    contDList[contID] = contItem;
  end
  return contList,contDList;
end
