using DelimitedFiles
using DataFrames
using CSV: read
include("NetworkData.jl")

"""
main parsing for all files and initialization of Power Network structure and Contingencies
in: pathes -- dictionary with pathes to files

TO DO
paralellize function calls. 
"""
function parser(pathes)
  # pathes to files
  rawFile = pathes[:raw]
  contFile = pathes[:contin]
  costsFile = pathes[:costs]
  inlFile = pathes[:inl]

  # create an empty Network structure and fill it with parsed data
  PN = PNetwork()

  # empty structures with DataFrames to store Power Network data and Contingencies data
  mainData = PNetworkData()
  costsData = CostsStruct()
  continData = ContingenciesStruct()
  # parsing
  raw_parser!(rawFile, mainData, PN); # fill the data into DataFrames and initialize PN with parameters based on case.raw
  response_parser!(inlFile, PN)
  costs_parser!(costsFile, costsData, PN)
  cont_parser!(contFile, continData)

  return PN, mainData, costsData, continData
end

"""
reads and returns file
in: path -- path to a file to read
out: file -- open file
"""
function file_reader(path)
  file = readdlm(path,',', skipblanks=false);
  return file
end


"""
pushes bunch of raw data into a dataframe
in: bunch -- data in "raw" format (delimited somehow)
    frame -- DataFrame to store the data
"""
function fillin_data!(bunch, frame)
    m,n = size(bunch)
    len = size(frame)[2]
    for i = 1:m
        push!(frame, bunch[i, 1:len])
    end
end

"""
custom function to store transformer data, since the data is strange -- one line of data splitted into 4 lines in the *.raw file.
in:  bunch -- data in "raw" format (delimited somehow)
     frame -- DataFrame to store the data
     type_of_data -- string parametr to specify the data is transformer data
"""
function fillin_data!(bunch, frame, type_of_data::String)
  m,n = size(bunch)
  frame_len = size(frame)[2]
  step = 4 # these 4 lines should be contacatinated as one to put into the frame
  lines_len = [21, 3, 17, 2] # length of each of 4 lines. 43 attributes of transformer in total
  for i = 1:step:m
    full_line = []
    for j = 0:step-1
      full_line = [full_line; bunch[i+j, 1:lines_len[j+1]]]; # concatination off 4 lines into 1
    end
    push!(frame, full_line)
  end
end


"""
reads *.raw data file, parses into bunches according to the data types (bus, load, branches and so on) and stores all the data into NetData structure inside corresponding DataFrames;
Initialize all parameters of PN
in: file_path -- path to .raw file
    NetData -- structure with empty DataFrames to store data
    PN -- structure representing Power Network
"""
function raw_parser!(file_path, NetData, PN)
    # rawData = readdlm(file_path,',', skipblanks=false);
    rawData = file_reader(file_path)

    n,m = size(rawData);
    # indeces of starting and ending lines for each bunch
    # bus
    busStartL = 4; # first three lines are headers
    busEndL = 0;
    # load
    loadStartL = 0;
    loadEndL = 0;
    # fixed shunt
    fixedShuntStartL = 0;
    fixedShuntEndL = 0;
    # generator
    genStartL = 0;
    genEndL = 0;
    # line
    branchStartL = 0;
    branchEndL = 0;
    # transfromer
    transformerStartL = 0;
    transformerEndL = 0;
    # switched shunt
    switchedShuntStartL = 0;
    switchedShuntEndL = 0;

    # parsing line by line to find indeces of each bunch
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
      end
      if occursin("END OF SWITCHED SHUNT DATA", string(rawData[i,1]))
          switchedShuntEndL = i - 1;
      end
    end

    # base MVA
    # NetData.sbase = rawData[1,2];
    PN.s = rawData[1,2];

    # bus
    fillin_data!(rawData[busStartL:busEndL,:], NetData.bus) # store into DataFrame
    Bus_init!(rawData[busStartL:busEndL,:], PN) # store into Array of Buses

    # load
    # println(loadStartL, loadEndL)
    fillin_data!(rawData[loadStartL:loadEndL,:], NetData.load)
    Load_init!(rawData[loadStartL:loadEndL, :], PN)

    # fixed shunt
    # current challenge does not have data for fixed shunt
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

end


"""
parser for costs. parse costs from file at costs_path into DataFrames of sctructure costsData used to initialize costs in PN.
steps: 1. Store all data into DataFrames
       2. Join corresponding DataFrames on described (in provided file) columns
       3. Using resulting new DataFrame initialize costs in PN

in: costs_path -- path to file with costs
    costsData -- structure with DataFrames
    PN -- Power Network structure


TO DO
modify if's to check for end of line instead of only "0 / " substring

"""
function costs_parser!(costs_path, costsData, PN)
  costs = file_reader(costs_path)

  n, m = size(costs)
  # STEP 1.  Fill data from file into DataFrames
  i = 1;
  while i <= n

    # Generator Dispatch Units
    if occursin("Generator Dispatch Units", string(costs[i,1]))
      i += 1
      len = size(costsData.genDispatch)[2]
      while !(occursin("0 / ", string(costs[i,1])))
        push!(costsData.genDispatch, costs[i, 1:len]) # store generator dispatch data
        i += 1
      end
      i -= 1
    end

    # Active Power Dispatch Tables
    if occursin("Active Power Dispatch Tables", string(costs[i,1]))
      i += 1
      len = size(costsData.activeDispatch)[2]
      while !(occursin("0 / ", string(costs[i,1])))
        push!(costsData.activeDispatch, costs[i, 1:len])
        i += 1
      end
      i -= 1
    end

    # Piece-wise Linear Cost Tables
    if occursin("Piece-wise Linear Cost Tables", string(costs[i,1]))
      i += 1
      len = size(costsData.linearCurve)[2] # numb. of columns
      while !(occursin("0 / ", string(costs[i,1])))
        line = costs[i,1:len-1] # first line (out of 4) with *LTBL, ‘LABEL’, *NPAIRS
        npairs = line[end] # indicates the number of lines with X-Y pairs
        XYbunch = costs[i+1:i+npairs, 1:2]
        lineToPush = [line[1], line[2], line[3], XYbunch] # *LTBL, ‘LABEL’, *NPAIRS, X-Y data
        push!(costsData.linearCurve, lineToPush)
        i += npairs+1
      end
      i -= 1
    end
    i += 1
  end

  # initialize costs in the network using above DataFrame
  costs_init!(costsData, PN)
end


"""
read and parse response factor α from .inl file directly into PN
in: resp_path -- path to .inl file
    PN -- Power Network structure
"""
function response_parser!(resp_path, PN)#, Resp)
  to_skip = 2 # there are to lines at the end "0" and "Q" which indicate the end of the file
  # Resp.data = read(resp_path, footerskip=to_skip, header=names(Resp.data))
  # Data = read(resp_path, footerskip=to_skip, header=names(Resp.data))
  open(resp_path) do file
    for L in eachline(file)
      if length(split(L)) == 1
        break
      end
      line = split(L) # split line to create generator index g = (i, id) and response factor
      i = parse(Int64, replace(line[1], ","=> ""))
      id = parse(Int64, replace(line[2], ","=> ""))
      alpha = parse(Float64, replace(line[6], ","=> ""))
      #
      PN.GeneratorList[PN.gen_ind_I[(i, id)]].alpha = alpha
    end
  end
end


"""
parser for contingencies. parses into structure contData
in: contin_path -- path to file .con
    contData -- strucure to store contingencies
"""
function cont_parser!(contin_path, contData)
  # gen_cont = [];
  # line_cont = [];

  open(contin_path) do file
    for L in eachline(file)
      # generator is off
      if occursin("REMOVE UNIT", L)
        push!(contData.genCont, split(L)[end])
      end
      # line is off
      if occursin("OPEN BRANCH", L)
        from = parse(Int, split(L)[end-5])
        to = parse(Int, split(L)[end-2])
        push!(contData.lineCont, (from, to))
      end
    end
  end
end


# """
# main parsing for all files and initialization of Power Network structure and Contingencies
# in: pathes -- dictionary with pathes to the files
# """
# function parser(pathes)
#   # create an empty Network structure and fill it with parsed data
#   PN = PNetwork()
#
#   # empty structures with DataFrames to store Power Network data and Contingencies data
#   mainData = PNetworkData()
#   costsData = CostsStruct()
#   continData = ContingenciesStruct()
#   # parsing
#   raw_parser!(rawFile, mainData, PN); # fill the data into DataFrames and initialize PN with parameters based on case.raw
#   response_parser!(inlFile, PN)
#   costs_parser!(costsFile, costsData, PN)
#   cont_parser!(contFile, continData)
#
#   return PN, mainData, costsData, continData
# end
