include("dProc.jl");
include("buildMod.jl");
include("dPrint.jl");

function MyJulia1(rawFile, genFile, contFile)

    println("RAW file: ",rawFile)
    println("gen file: ",genFile)
    println("con file: ",contFile)

    println("Data generation")

    tic()

    BLGS, Br, contingency, genSeg, baseMVA = dProc(rawFile, genFile, contFile);
    toc()

    println("Model generation")

    tic()
    PSCOPF = buildMod(BLGS, Br, contingency, genSeg);
    toc()

    println("Model solve")
    tic()
    solve(PSCOPF);
    toc()

    println("Solution writing")
    tic()
    dPrint(PSCOPF, genSeg, BLGS, Br, contingency, baseMVA);
    toc()
end
