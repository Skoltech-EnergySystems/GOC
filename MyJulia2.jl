include("def.jl");
include("dProc.jl");
include("buildMod.jl");
include("dPrint.jl");

function MyJulia2(rawFile, genFile, contFile)

    println("RAW file: ",rawFile)
    println("gen file: ",genFile)
    println("con file: ",contFile)

    println("Data generation")

    tic()
    baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg = preProc(rawFile);
    busList,busDList,genList,genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
    brList,brListSingle,brDList = branchProc(baseMVA,branchSeg,transformerSeg);
    contingency,contList,contDList = contProc(contFile);
    fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
    uData = uncertainData(contList,contDList);

    toc()

    println("Model generation")

    tic()
    PSCOPF,NGen,NBus,NBr,NK,aL = buildMod(fData,uData, contDList, contingency);
    toc()

    println("Model solve")
    tic()
    solve(PSCOPF);
    toc()

    println("Solution writing")
    tic()
    dPrint(PSCOPF,NGen,NBus,NBr,NK,aL,genSeg,busSeg,brList,baseMVA);
    toc()
end
