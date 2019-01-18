include("DataParsing.jl")

rawFile = "./data/case.raw"
baseMVA, busSeg, loadSeg, shuntSeg, genSeg, branchSeg, transformerSeg = preProc(rawFile);




function contProc2(contFile)
  # path = "./data/" * contFile
  gen_cont = [];
  line_cont = [()];

  open(contFile) do file
    for L in eachline(file)
      if contains(L, "REMOVE UNIT")
        push!(gen_cont, split(L)[end])
      end
      if contains(L, "OPEN BRANCH")
        from = parse(split(L)[end-5])
        to = parse(split(L)[end-2])
        push!(line_cont, (from, to))
      end
    end
  end

  return gen_cont, line_cont
end



contFile = "./data/case.con";
a, b = contProc2(contFile)
