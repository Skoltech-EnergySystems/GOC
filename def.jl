# definitions of the data type

type busData
  ID :: Int64
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
  ID :: Any
  Name :: Any
  Loc :: Int64

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
  ID :: Any
  revID :: Any

  r :: Float64
  x :: Float64
  g :: Float64
  b :: Float64
  bc :: Float64

  t :: Float64
  zeroImpe :: Bool
  tau :: Float64
  tauprime :: Float64
  thetatr :: Float64
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
