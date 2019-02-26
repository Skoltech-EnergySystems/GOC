using DelimitedFiles
using DataFrames
using CSV: read
include("NetworkData.jl")

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
function raw_parser!(file_path, NetData, PN)
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
      if occursin("END OF BUS DATA", string(rawData[i,1]))
        busEndL = i - 1;
        loadStartL = i + 1;
      end
      if occursin("END OF LOAD DATA", string(rawData[i,1]))
        loadEndL = i - 1;
        fixedShuntStartL = i + 1;
      end
      if occursin("END OF FIXED SHUNT DATA", string(rawData[i,1]))
        fixedShuntEndL = i - 1;
        genStartL = i + 1;
      end
      if occursin("END OF GENERATOR DATA", string(rawData[i,1]))
        genEndL = i - 1;
        branchStartL = i + 1;
      end
      if occursin("END OF BRANCH DATA", string(rawData[i,1]))
        branchEndL = i - 1;
        transformerStartL = i + 1;
      end
      if occursin("END OF TRANSFORMER DATA", string(rawData[i,1]))
        transformerEndL = i - 1;
      end
      if occursin("END OF FACTS CONTROL DEVICE DATA", string(rawData[i,1]))
        switchedShuntStartL = i + 1
        # println(switchedShuntStartL)
      end
      if occursin("END OF SWITCHED SHUNT DATA", string(rawData[i,1]))
          switchedShuntEndL = i - 1;
          # println(switchedShuntEndL)
      end
    end

    # base MVA
    # NetData.sbase = rawData[1,2];
    PN.s = rawData[1,2];

    # bus
    fillin_data!(rawData[busStartL:busEndL,:], NetData.bus) # store into DataFrame
    Bus_init!(rawData[busStartL:busEndL,:], PN) # store into Array of Buses

    # load
    println(loadStartL, loadEndL)
    fillin_data!(rawData[loadStartL:loadEndL,:], NetData.load)
    Load_init!(rawData[loadStartL:loadEndL, :], PN)

    # fixed shunt
    fillin_data!(rawData[fixedShuntStartL:fixedShuntEndL,:], NetData.fixedBusShunt)
    if fixedShuntEndL > fixedShuntStartL + 1
      fShunt_init!(rawData[fixedShuntStartL:fixedShuntEndL,:], PN)
    end

    # generator
    fillin_data!(rawData[genStartL:genEndL,:], NetData.generator)
    Generator_init!(rawData[genStartL:genEndL,:], PN)

    # branch LINE
    fillin_data!(rawData[branchStartL:branchEndL,:], NetData.line)
    Line_init!(rawData[branchStartL:branchEndL,:], PN)

    # transformer
    fillin_data!(rawData[transformerStartL:transformerEndL,:], NetData.transformer, "transformer")
    Transformer_init!(rawData[transformerStartL:transformerEndL,:], PN)

    # switched shunt
    fillin_data!(rawData[switchedShuntStartL:switchedShuntEndL, :], NetData.switchedShunt)
    sShunt_init!(rawData[switchedShuntStartL:switchedShuntEndL, :], PN)

    # println("raw data at $file_path is parsed")
end


"""
parser for costs

TO DO
modify if's to check for end of line instead of only "0 / " substring

"""
function costs_parser!(costs_path, costsData, PN)
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
    if occursin("Generator Dispatch Units", string(costs[i,1]))
      # println("\t=====================\n", "\t $i | inside first if\n")
      # println("\t with string | $(costs[i,1]) |")
      # println("\t=====================")
      # s[1] += 1
      i += 1
      len = size(costsData.genDispatch)[2]
      while !(occursin("0 / ", string(costs[i,1])))
        # println("\t\t $i | looping in while in 1 if with string | $(costs[i,1])")
        push!(costsData.genDispatch, costs[i, 1:len])
        i += 1
      end
      i -= 1
    end

    # Active Power Dispatch Tables
    if occursin("Active Power Dispatch Tables", string(costs[i,1]))
      # println("\t=====================\n", "\t $i | inside second if\n")
      # println("\t  with string | $(costs[i,1]) |")
      # println("\t=====================")
      # s[2] += 1
      i += 1
      len = size(costsData.activeDispatch)[2]
      while !(occursin("0 / ", string(costs[i,1])))
        # println("\t\t $i | looping in while in 2 if with string | $(costs[i,1])")
        # println(string(costs[i+1,:]))
        push!(costsData.activeDispatch, costs[i, 1:len])
        i += 1
      end
      i -= 1
    end

    # Piece-wise Linear Cost Tables
    if occursin("Piece-wise Linear Cost Tables", string(costs[i,1]))
      # println("\t=====================\n", "\t $i | inside third if\n")
      # println("\t with string | $(costs[i,1]) |")
      # println("\t=====================")
      # s[3] += 1
      i += 1
      len = size(costsData.linearCurve)[2]
      while !(occursin("0 / ", string(costs[i,1])))
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

  # initialize costs in the network
  costs_init!(costsData, PN)
end


# RESPONSE in case.inl
function response_parser!(resp_path, PN)#, Resp)
  to_skip = 2 # there are to lines at the end "0" and "Q" which indicate the end of the file
  # Resp.data = read(resp_path, footerskip=to_skip, header=names(Resp.data))
  # Data = read(resp_path, footerskip=to_skip, header=names(Resp.data))
  open(resp_path) do file
    for L in eachline(file)
      if length(split(L)) == 1
        break
      end
      line = split(L)
      i = parse(Int64, replace(line[1], ","=> ""))
      id = parse(Int64, replace(line[2], ","=> ""))
      alpha = parse(Float64, replace(line[6], ","=> ""))
      #
      PN.GeneratorList[PN.gen_ind_I[(i, id)]].alpha = alpha
    end
  end
end


"""
parser for contingencies
"""
function cont_parser!(contin_path, contData)
  # gen_cont = [];
  # line_cont = [];

  open(contin_path) do file
    for L in eachline(file)
      if occursin("REMOVE UNIT", L)
        push!(contData.genCont, split(L)[end])
      end
      if occursin("OPEN BRANCH", L)
        from = parse(Int, split(L)[end-5])
        to = parse(Int, split(L)[end-2])
        # line = (from, to)
        push!(contData.lineCont, (from, to))
      end
    end
  end
end


  # println("contingencies at $contin_path are parsed ")
