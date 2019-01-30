using DataFrames
using CSV: read

"""
TO DO
input: all files to parse
output: structure with all neede data stored in DataFrames
in particular should contain calls of 4 parser for raw, contingencies, costs and response files
"""
function main_parser()

end

"""
reads and returns file
"""
function file_reader(path)
  file = readdlm(path,',', skipblanks=false);
  return file
end


"""
pushes bunch of raw data into a dataframe
"""
function fillin_data!(bunch, frame)
    m,n = size(bunch)
    len = size(frame)[2]
    for i = 1:m
        push!(frame, bunch[i, 1:len])
    end
end

"""
custom function to store transformer data, since the data is strange -- one line of data splitted into 4 lines in the *.raw file. it may possibly change further.
"""
function fillin_data!(bunch, frame, type_of_data::String)
  m,n = size(bunch)
  frame_len = size(frame)[2]
  step = 4 # these 4 lines should be contacatinated as one to put into the frame
  lines_len = [21, 3, 17, 2] # length of each of 4 lines. 43 attributes of transformer in total
  for i = 1:step:m
    full_line = []
    for j = 0:step-1
      full_line = [full_line; bunch[i+j, 1:lines_len[j+1]]];
    end
    push!(frame, full_line)
  end
end


"""
reads *.raw data file, parses into bunches according to the data types (bus, load, branches and so on) and stores all the data into NetData structure
"""
function raw_parser!(file_path, NetData)
    # rawData = readdlm(file_path,',', skipblanks=false);
    rawData = file_reader(file_path)

    n,m = size(rawData);
    busStartL = 4; # first three lines are headers
    # indeces of starting and ending lines for each bunch
    busEndL = 0;
    loadStartL = 0;
    loadEndL = 0;
    fixedShuntStartL = 0;
    fixedShuntEndL = 0;
    genStartL = 0;
    genEndL = 0;
    branchStartL = 0;
    branchEndL = 0;
    transformerStartL = 0;
    transformerEndL = 0;
    switchedShuntStartL = 0;
    switchedShuntEndL = 0;

    # parsing by line to find indeces of each bunch
    for i in 1:n
      if contains(string(rawData[i,1]),"END OF BUS DATA")
        busEndL = i - 1;
        loadStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF LOAD DATA")
        loadEndL = i - 1;
        fixedShuntStartL = i + 1;
      end
      if contains(string(rawData[i,1]),"END OF FIXED SHUNT DATA")
        fixedShuntEndL = i - 1;
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
      if contains(string(rawData[i,1]), "END OF FACTS CONTROL DEVICE DATA")
        switchedShuntStartL = i + 1
        # println(switchedShuntStartL)
      end
      if contains(string(rawData[i,1]), "END OF SWITCHED SHUNT DATA")
          switchedShuntEndL = i - 1;
          # println(switchedShuntEndL)
      end
    end

    """
    the next block can be changed to a function call which will store the data directly into PNetwork with all neede computations
      the same is true for the rest of parsers
    """
    # base MVA
    NetData.sbase = rawData[1,2];
    # bus
    fillin_data!(rawData[busStartL:busEndL,:], NetData.bus)
    # load
    fillin_data!(rawData[loadStartL:loadEndL,:], NetData.load)
    # fixed shunt
    fillin_data!(rawData[fixedShuntStartL:fixedShuntEndL,:], NetData.fixedBusShunt)
    # generator
    fillin_data!(rawData[genStartL:genEndL,:], NetData.generator)
    # branch
    fillin_data!(rawData[branchStartL:branchEndL,:], NetData.branch)
    # transformer
    fillin_data!(rawData[transformerStartL:transformerEndL,:], NetData.transformer, "transformer")
    # switched shunt
    fillin_data!(rawData[switchedShuntStartL:switchedShuntEndL, :], NetData.switchedShunt)

    println("raw data at $file_path is parsed")
end


"""
parser for contingencies
"""

function cont_parser!(contin_path, contData)
  # gen_cont = [];
  # line_cont = [];

  open(contin_path) do file
    for L in eachline(file)
      if contains(L, "REMOVE UNIT")
        push!(contData.genCont, split(L)[end])
      end
      if contains(L, "OPEN BRANCH")
        from = parse(split(L)[end-5])
        to = parse(split(L)[end-2])
        # line = (from, to)
        push!(contData.lineCont, (from, to))
      end
    end
  end
  println("contingencies at $contin_path are parsed ")
  # return gen_cont, line_cont
end


"""
parser for costs

TO DO
modify if's to check for end of line instead of only "0 / " substring

"""
function costs_parser!(costs_path, costsData)
  # costs = readdlm(costsFile,',', skipblanks=false);
  costs = file_reader(costs_path)

  n, m = size(costs)
  # s = [0 0 0];
  # println("starting loop")
  i = 1;
  while i <= n
    # i += 1
    # println("i=$i | ", costs[i,:])
    # Generator Dispatch Units
    if contains(string(costs[i,1]), "Generator Dispatch Units")
      # println("\t=====================\n", "\t $i | inside first if\n")
      # println("\t with string | $(costs[i,1]) |")
      # println("\t=====================")
      # s[1] += 1
      i += 1
      len = size(costsData.genDispatch)[2]
      while !(contains(string(costs[i,1]), "0 / "))
        # println("\t\t $i | looping in while in 1 if with string | $(costs[i,1])")
        push!(costsData.genDispatch, costs[i, 1:len])
        i += 1
      end
      i -= 1
    end

    # Active Power Dispatch Tables
    if contains(string(costs[i,1]), "Active Power Dispatch Tables")
      # println("\t=====================\n", "\t $i | inside second if\n")
      # println("\t  with string | $(costs[i,1]) |")
      # println("\t=====================")
      # s[2] += 1
      i += 1
      len = size(costsData.activeDispatch)[2]
      while !(contains(string(costs[i,1]), "0 / "))
        # println("\t\t $i | looping in while in 2 if with string | $(costs[i,1])")
        # println(string(costs[i+1,:]))
        push!(costsData.activeDispatch, costs[i, 1:len])
        i += 1
      end
      i -= 1
    end

    # Piece-wise Linear Cost Tables
    if contains(string(costs[i,1]), "Piece-wise Linear Cost Tables")
      # println("\t=====================\n", "\t $i | inside third if\n")
      # println("\t with string | $(costs[i,1]) |")
      # println("\t=====================")
      # s[3] += 1
      i += 1
      len = size(costsData.linearCurve)[2]
      while !(contains(string(costs[i,1]), "0 / "))
        # println("\t\t $i | looping in while in 3 if with string | $(costs[i,1])")
        line = costs[i,1:len-1] # title line with *LTBL, ‘LABEL’, *NPAIRS
        npairs = line[end] # indicates the number of lines with X-Y pairs
        XYbunch = costs[i+1:i+npairs, 1:2]
        lineToPush = [line[1], line[2], line[3], XYbunch]
        push!(costsData.linearCurve, lineToPush)
        i += npairs+1
      end
      i -= 1
    end
    i += 1
  end
  # println("end of loop")
  println("costs file at $costs_path is parsed")
  # return k, A, s
end


"""
it is broken now!!!
"""
# RESPONSE in case.inl
function response_parser!(resp_path, Resp)
  to_skip = 2 # there are to lines at the end "0" and "Q" which indicate the end of the file
  Resp.data = read(resp_path, footerskip=to_skip, header=names(Resp.data))
end


####################### HELPERS-TESTERS

# function analyze_seq(seq)
#   m, n = size(seq)
#   println(size(seq))
#   for i = 1:5
#     println(seq[i, 1:13][1:end])
#   end
# end
#
# function analyze(file_path)
#   rawData = readdlm(file_path,',', skipblanks=false);
#   n,m = size(rawData);
#   switchedShuntStartL = 0;
#
#   for i  = 1:n
#     if contains(string(rawData[i,1]), "END OF FACTS CONTROL DEVICE DATA")
#       switchedShuntStartL = i + 1
#       println(switchedShuntStartL)
#     end
#   end
# end
