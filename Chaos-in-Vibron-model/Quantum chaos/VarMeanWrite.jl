using Plots
using PyCall
using Printf
using LinearAlgebra
using StatsBase
using Statistics


#
#
#   writes var and mean value of OTOC and spectrum
#
#

include("SpectrumFunctions.jl")

function OTOCevaluation(W,V,t::Float64,spectrum::Array{Float64,1})

    Dim = length(spectrum)
    W_t = Diagonal(exp.( 1im *t*spectrum)) * W * Diagonal(exp.( -1im *t*spectrum))

    WV = W_t * V
    VW = V * W_t

    OTOC = 2/Dim * ((WV - VW)*VW - (WV - VW)*WV)

    OTOCvalues = Array{Float64,1}(undef,Dim)
    for j in 1:Dim
        OTOCvalues[j] = OTOC[j,:]'*OTOC[j,:]
    end

    return OTOCvalues
end


function OTOCVarMean(W, V, Hamiltonian; len::Int64 = 1500, Min = 1e7, Max = 1e8)
    

    EigSystem = eigen(Hamiltonian)
    spectrum = EigSystem.values
    vectors = EigSystem.vectors
    
    Dim = length(spectrum)
    
    S = vectors
    W_matrix = inv(S)*W*S #n ve vlastni bazi hamiltonianu
    V_matrix = inv(S)*V*S
    

    OTOCarrayEnergyScale = Array{Float64,2}(undef,len,Dim)
    MeanEnergyScale = Array{Float64,1}(undef,Dim)
    VarEnergyScale = Array{Float64,1}(undef,Dim)


    Points = sort(rand(len))
    Interval = (Max - Min)

    Threads.@threads for i in 1:len
            println("$i" )
             
            OTOCvalues = OTOCevaluation(W_matrix,V_matrix,Min + Points[i]*Interval,spectrum)

            for j in 1:Dim
                OTOCarrayEnergyScale[i,j] = OTOCvalues[j]

            end
            
        end

        for j in 1:Dim
            MeanEnergyScale[j] = mean(OTOCarrayEnergyScale[:,j])
            VarEnergyScale[j] = sqrt(var(OTOCarrayEnergyScale[:,j]))
       
        end

    return MeanEnergyScale, VarEnergyScale, spectrum
end


function WriteOTOCFast!(N::Int64,ξ::Float64,ϵ::Float64;path::String = "")

    
    Qp = Qp_Nnl(N)
    Qm = Qm_Nnl(N)
    
    Pp = Pp_Nnl(N)
    Pm = Pm_Nnl(N)
    
    qx = -sqrt(2)*(Qp - Qm)
    px = -sqrt(2)*(Pp - Pm)
    OTOCname = "[qx(t),px]"
    
    H = Hamiltonian_Nnl(ξ,ϵ,N)
  

    MeanEnergy, VarEnergy, spectrum = OTOCVarMean(qx,px,H)
    

    pars = @sprintf "_N=%1.0i,chi=%1.2f,eps=%1.2f" N ξ ϵ

    open(path*"Spectrum"*OTOCname*pars*".txt", "a") do file
        println(file,spectrum)
    end

    open(path*"Mean"*OTOCname*pars*".txt", "a") do file
        println(file,MeanEnergy)
    end


    open(path*"Var"*OTOCname*pars*".txt", "a") do file
        println(file,VarEnergy)
    end



end
    


WriteOTOCFast!(100,0.5,0.1)





