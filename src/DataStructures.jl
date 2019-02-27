using DataFrames

#=
There are structures to store data from case files into DataFrames and their conctructors
All these functions are to create empty DataFrames with column names described in the instructions file
=#

# structure to store Network Data
mutable struct PNetworkData
    # sbase::Float16
    bus::DataFrame
    load::DataFrame
    fixedBusShunt::DataFrame
    generator::DataFrame
    line::DataFrame
    transformer::DataFrame
    switchedShunt::DataFrame

    PNetworkData() = constr_PNetworkData(new())

end

function constr_PNetworkData(PND)
    # BUS
    PND.bus = DataFrame(I=Int[], NAME=String[], BASKV=Float16[], IDE=Float16[], AREA=Int[],           ZONE=Float16[], OWNER=Float16[], VM=Float16[], VA=Float16[], NVHI=Float16[], NVLO=Float16[], EVHI=Float16[], EVLO=Float16[])
    #
    # PND.bus = busData

    # LOAD
    PND.load = DataFrame(I=Int[], ID=String[], STATUS=Bool[], AREA=Float16[], ZONE=Float16[], PL=Float16[], QL=Float16[], IP=Float16[], IQ=Float16[], YP=Float16[], YQ=Float16[], OWNER=Float16[], SCALE=Float16[], INTRPT=Float16[])
    #
    # PND.load = loadData

    # Fixed Bus ShunT
    PND.fixedBusShunt = DataFrame(I=Int[], ID=String[], STATUS=Bool[], GL=Float16[], BL=Float16[])
    #
    # PND.fixedBusShunt = fixedBusShunt

    # GENERATOR
    PND.generator = DataFrame(I=Int[], ID=String[], PG=Float16[], QG=Float16[], QT=Float16[], QB=Float16[], VS=Float16[], IREG=Float16[], MBASE=Float16[], ZR=Float16[], ZX=Float16[], RT=Float16[], XT=Float16[], GTAP=Float16[], STAT=Bool[], RMPCT=Float16[],
    PT=Float16[], PB=Float16[], O1=Float16[], F1=Float16[], O2=Float16[], F2=Float16[], O3=Float16[], F3=Float16[], O4=Float16[], F4=Float16[], WMOD=Float16[], WPF=Float16[])
    #
    # PND.generator = generator

    # LINE
    PND.line = DataFrame(I=Int[], J=Int[], CKT=String[], R=Float16[], X=Float16[], B=Float16[], RATEA=Float16[], RATEB=Float16[], RATEC=Float16[], GI=Float16[], BI=Float16[], GJ=Float16[], BJ=Float16[], ST=Bool[], MET=Float16[], LEN=Float16[], O1=Float16[], F1=Float16[], O2=Float16[], F2=Float16[], O3=Float16[], F3=Float16[], O4=Float16[], F4=Float16[])
    #
    # PND.line = line

    # TRANSFORMER
    PND.transformer = DataFrame(I=Int[], J=Int[], K=Float16[], CKT=String[], CW=Float16[], CZ=Float16[], CM=Float16[], MAG1=Float16[], MAG2=Float16[], NMETR=Float16[], NAME=String[], STAT=Bool[], O1=Float16[], F1=Float16[], O2=Float16[], F2=Float16[], O3=Float16[], F3=Float16[], O4=Float16[], F4=Float16[], VECGRP=String[],
    R12=Float16[], X12=Float16[], SBASE12=Float16[],
    WINDV1=Float16[], NOMV1=Float16[], ANG1=Float16[], RATA1=Float16[], RATB1=Float16[], RATC1=Float16[], COD1=Float16[], CONT1=Float16[], RMA1=Float16[], RMI1=Float16[], VMA1=Float16[], VMI1=Float16[], NTP1=Float16[], TAB1=Float16[], CR1=Float16[], CX1=Float16[], CNXA1=Float16[],
    WINDV2=Float16[], NOMV2=Float16[])
    #
    # PND.transformer = transformer

    # Switched Shunt
    PND.switchedShunt = DataFrame(I=Int[], MODSW=Float16[], ADJM=Float16[], STAT=Bool[], VSWHI=Float16[], VSWLO=Float16[], SWREM=Float16[], RMPCT=Float16[], RMIDNT=String[], BINIT=Float16[], N1=Int[], B1=Float16[], N2=Int[], B2=Float16[], N3=Int[], B3=Float16[], N4=Int[], B4=Float16[], N5=Int[], B5=Float16[], N6=Int[], B6=Float16[], N7=Int[], B7=Float16[], N8=Int[], B8=Float16[])
    #
    # PND.switchedShunt = switchedShunt

    return PND
end

###### CONTINGENCIES
mutable struct ContingenciesStruct
    genCont::Array
    lineCont::Array

    ContingenciesStruct() = constr_ContingenciesStruct(new())
end
"""
ContingenciesStruct constructor
"""
function constr_ContingenciesStruct(C)
    C.genCont = []
    C.lineCont = []
    return C
end


######## COSTS
mutable struct CostsStruct
    genDispatch::DataFrame
    activeDispatch::DataFrame
    linearCurve::DataFrame

    CostsStruct() = initialize_costs_struct(new())
end
"""
CostsStruct constructor
"""
function initialize_costs_struct(C)
    C.genDispatch = DataFrame(BUS=Int[], GENID=Int[], DISP=Float16[], DSPTBL=Int[])

    #
    C.activeDispatch = DataFrame(TBL=Int[], PMAX=Float16[], PMIN=Float16[], FUELCOST=Float16[], CTYP=Float16[], STATUS=Float16[], CTBL=Int[])
    #
    C.linearCurve = DataFrame(LTBL=Int[], LABEL=String[], NPAIRS=Int[], XYmatrix=Array{Float16, 2}[])

    # genCosts = CostsStruct(genDispatch, activeDispatch, linearCurve)

    # rename for further join on these columns
    rename!(C.genDispatch, Dict(:DSPTBL => :TBL))
    rename!(C.linearCurve, Dict(:LTBL => :CTBL))

    return C
end

# RESPONSE FROM case.inl
mutable struct ResponseStruct
    data::DataFrame
end
"""
ResponseStruct constructor
"""
function initialize_response_struct()
    data = DataFrame(I=Int16[], ID=Int[], H=Float16[], PMAX=Float16[], PMIN=Float16[], R=Float16[], D=Float16[])
    Resp = ResponseStruct(data)
    return Resp
end
