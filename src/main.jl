include("DataParsing.jl")

rawFile = "./data/case.raw"
baseMVA, busSeg, loadSeg, shuntSeg, genSeg, branchSeg, transformerSeg = preProc(rawFile);
# process buses data
busList, busDList, genList, genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
# process branch data
brList, brListSingle, brDList = branchProc(baseMVA,branchSeg,transformerSeg);
# process contingencies
contingency, contList, contDList = contProc(contFile);



function contProc2(contFile)
  gen_cont = [];
  line_cont = [];

  open(contFile) do file
    for L in eachline(file)
      if contains(L, "REMOVE UNIT")
        push!(gen_cont, split(L)[end])
      end
      if contains(L, "OPEN BRANCH")
        from = parse(split(L)[end-5])
        to = parse(split(L)[end-2])
        # line = (from, to)
        push!(line_cont, (from, to))
      end
    end
  end

  return gen_cont, line_cont
end



# contFile = "./data/case.con";
# a, b = contProc2(contFile);



# contFile = "IEEE14-1_contingency.csv";
# contingency, contList, contDList = contProc(contFile);
