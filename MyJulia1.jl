include("def.jl");
include("dProc.jl");
include("buildMod.jl");
include("dPrint.jl");
#include("dPrint_pl.jl");

function MyJulia1(rawFile, genFile, contFile)

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
    PSCOPF,NGen,NBus,NBr,NBrc,NK,aL,Ct,Ctc,Cf,Cfc,Yf,Yfc,Yt,Ytc,br_int,br_intc = buildMod(fData,uData, contDList, contingency);
    toc()

    println("Model solve")
    tic()
    solve(PSCOPF);
    toc()

    println("Solution writing")
    tic()
    dPrint(PSCOPF,NGen,NBus,NBr,NBrc,NK,aL,Ct,Ctc,Cf,Cfc,Yf,Yfc,Yt,Ytc,br_int,br_intc,genSeg,busSeg,brList,brDList,baseMVA);
    #dPrint_pl(PSCOPF,NGen,NBus,NBr,NK,aL,Ct,Cf,Yf,Yt,genSeg,busSeg,brList,brDList,baseMVA);
    toc()
end
