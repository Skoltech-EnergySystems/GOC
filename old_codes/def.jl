# definitions of the data type

type busData
    ID_ext :: Int64
    ID_int :: Int64
    Name :: String
    gen :: Array{Any,1}

    Vmax :: Float64
    Vmin :: Float64
    Pd :: Float64
    Qd :: Float64

    gsh :: Float64
    bsh :: Float64
 end

type genData
    ID_ext :: Any
    Name :: Any
    Loc :: Int64
    ID_int :: Any

    cn :: Array{Int64,1}
    cParams :: Dict{Any,Any}

    Pmax :: Float64
    Pmin :: Float64
    Qmax :: Float64
    Qmin :: Float64

    alpha :: Float64
end

type branchData
    From :: Int64
    To :: Int64
    CKT :: Any
    ID_ext :: Any
    #revID :: Any

    r :: Float64
    x :: Float64
    bc :: Float64

    t :: Float64
    tau :: Float64
 end

type fixedData
  baseMVA :: Float64

  busList :: Array{Any,1}
  busDList :: Dict{Any,Any}

  genList :: Array{Any,1}
  genDList :: Dict{Any,Any}

  brList :: Array{Any,1}
  brListSingle :: Array{Any,1}
  brDList :: Dict{Any,Any}
end

type contingencyData
  ID :: Int64
  Type :: String
  Loc :: Array{Any,1}
end

type uncertainData
  contList :: Array{Int64,1}
  contDList :: Dict{Any,Any}
end
