include("DataParsing.jl")

rawFile = "./data/case.raw";
contFile = "./data/case.con";
genFile = "./data/case.rop";

# pre-processing; returns rawdata bunches for generation, branches, contingencies, etc.
baseMVA, busSeg, loadSeg, shuntSeg, genSeg, branchSeg, transformerSeg = preProc(rawFile);
# process buses data
busList, busDList, genList, genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
# process branch data
brList, brListSingle, brDList = branchProc(baseMVA,branchSeg,transformerSeg);
# process contingencies
# contingency, contList, contDList = contProc(contFile);

# intialize Power System parameters
fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
# initialize contingencies
uData = uncertainData(contList,contDList);
