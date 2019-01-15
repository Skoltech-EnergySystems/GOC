# preprocessing: read in the data and divide it into segments
function preProc(rawFile)

    rawData = readdlm(rawFile,',', skipblanks=false);
    n,m = size(rawData);
    busStartL = 4; # first three lines are headers
    busEndL = 0;
    loadStartL = 0;
    loadEndL = 0;
    shuntStartL = 0;
    shuntEndL = 0;
    genStartL = 0;
    genEndL = 0;
    branchStartL = 0;
    branchEndL = 0;
    transformerStartL = 0;
    transformerEndL = 0;
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

  busSeg_int = collect(1:bn);
  ### BUS DATA
  busList = [];
  busDList = Dict();
  busNo = 0;
  for i in 1:bn
    # build the bus item
    busID_ext = busSeg[i,1];
    busID_int = busSeg_int[i,1];
    push!(busList,busID_ext);
    busName = string(busSeg[i,2]);
    busVmax = busSeg[i,10];
    busVmin = busSeg[i,11];
    busItem = busData(busID_ext,busID_int,busName,[],busVmax,busVmin,0,0,0,0);
    busDList[busID_ext] = busItem;
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
    busID_ext = genSeg[i,1];
    NGen = size(genSeg,1);
    genID_int = (1:NGen)[i];
    genLoc = busID_ext;
    genName = genSeg[i,2];
    genID = (genLoc,genName,genID_int);
    push!(genList,genID);
    push!(busDList[busID_ext].gen,genID);

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

    genItem = genData(busID_ext,genName,genLoc,genID_int,genCn,gencParams,genPmax,genPmin,genQmax,genQmin,genalpha);
    genDList[genID] = genItem;
  end
  return busList,busDList,genList,genDList;
  end

  function branchProc(baseMVA,branchSeg,transformerSeg)
  brn,brm = size(branchSeg);
  brList = [];
  brListSingle = [];
  brDList = Dict();
  for i in 1:brn
    brFrom = branchSeg[i,1];
    brTo = branchSeg[i,2];
    brCKT = branchSeg[i,3];
    branchID_ext = (brFrom,brTo,brCKT);

    brr = branchSeg[i,4];
    brx = branchSeg[i,5];
    brbc = branchSeg[i,6];
    brt = branchSeg[i,7]/baseMVA;

    brtau = 1;
    brItem = branchData(brFrom,brTo,brCKT,branchID_ext,brr,brx,brbc,brt,brtau);

    push!(brList,branchID_ext);
    push!(brListSingle,branchID_ext);
    brDList[branchID_ext] = brItem;
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
    trID_ext = (trFrom,trTo,trName);

    trr = transformerSeg[lineNo+1,1];
    trx = transformerSeg[lineNo+1,2];
    trbc = 0;
    trt = transformerSeg[lineNo+2,4]/baseMVA;

    trtau = transformerSeg[lineNo+2,1]/transformerSeg[lineNo+3,1];
    trItem = branchData(trFrom,trTo,trName,trID_ext,trr,trx,trbc,trt,trtau);

    push!(brList,trID_ext);
    push!(brListSingle,trID_ext);
    brDList[trID_ext] = trItem;
    lineNo = lineNoNew;
  end
  return brList,brListSingle,brDList;
  end

function contProc(contFile)
contingency = CSV.read(contFile);
contData,contTitle = readdlm(contFile,',',header = true);
contn,contm = size(contData);
contList = [];
contDList = Dict();
for i in 1:contn
  contID = contData[i,1];
  push!(contList,contID);
  contType = contData[i,2];
  if (contType == "B")|(contType == "T")
    contLoc = [(contData[i,3],contData[i,4],contData[i,5])];
  else
    contLoc = [contData[i,3],contData[i,4]];
  end
  contItem = contingencyData(contID,contType,contLoc);
  contDList[contID] = contItem;
end
return contingency,contList,contDList;
end
