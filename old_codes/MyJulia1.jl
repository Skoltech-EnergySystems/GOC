include("def.jl");
include("dProc.jl");
include("buildMod.jl");
include("dPrint.jl");

function MyJulia1(rawFile, genFile, contFile)

    println("RAW file: ",rawFile)
    println("gen file: ",genFile)
    println("con file: ",contFile)

    println("Data generation")

    tic()
    # pre-processing; returns rawdata bunches for generation, branches, contingencies, etc.
    baseMVA, busSeg, loadSeg, shuntSeg, genSeg, branchSeg, transformerSeg = preProc(rawFile);
    # process buses data
    busList, busDList, genList, genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
    # process branch data
    brList, brListSingle, brDList = branchProc(baseMVA,branchSeg,transformerSeg);
    # process contingencies
    contingency, contList, contDList = contProc(contFile);

    # intialize Power System parameters
    fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
    # initialize contingencies
    uData = uncertainData(contList,contDList);

    toc()

    println("Model generation")

    tic()
    PSCOPF,NGen,NBus,NBr,NK,aL,Ct,Cf,Yf,Yt = buildMod(fData, uData, contDList, contingency);
    toc()

    println("Model solve")
    tic()
    solve(PSCOPF);
    toc()

    println("Solution writing")
    tic()
    dPrint(PSCOPF,NGen,NBus,NBr,NK,aL,Ct,Cf,Yf,Yt,genSeg,busSeg,brList,brDList,baseMVA);
    toc()
end
