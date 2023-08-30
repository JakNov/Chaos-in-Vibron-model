using Plots
using PyCall
using Printf
using LinearAlgebra
using Statistics
using LsqFit
using SpecialFunctions
using StatsBase
using Polynomials


include("SpectrumFunctions.jl")


function GOE(n::Int)
    A = randn(n, n)
    normalization = 1 / √(2n)
    return Symmetric((A + A') / 2 * normalization)
end

function OTOCevaluation(W,V,t::Float64,spectrum::Array{Float64,1})

    Dim = length(spectrum)
    #println("diagonal multiplication time")
    W_t = Diagonal(exp.( 1im *t*spectrum)) * W * Diagonal(exp.( -1im *t*spectrum))

    WV = W_t * V
    VW = V * W_t

    OTOC = #=1/Dim *=# (WV-VW)

    OTOCvalues = Array{Float64,1}(undef,Dim)
    for j in 1:Dim
        OTOCvalues[j] = OTOC[j,:]'*OTOC[j,:]
    end

    return OTOCvalues
end


function OTOCgeneralRegime(W, V, Hamiltonian; len = 10000,tmax = 1e9,tmin = 1e7)
    

    EigSystem = eigen(Hamiltonian)
    spectrum = EigSystem.values
    vectors = EigSystem.vectors
    
    Dim = length(spectrum)
    
    S = vectors
    W_matrix = inv(S)*W*S #n ve vlastni bazi hamiltonianu
    V_matrix = inv(S)*V*S
    

    OTOCarrayEnergyScale = Array{Float64,2}(undef,len,Dim)
    times = Array{Float64,1}(undef,len)

    ##OTOCarrayParticleScale = Array{Float64,2}(undef,len,Dim)


    Threads.@threads for i in 1:len
            println("$i" )
            time = tmin + (tmax - tmin)*(i-1)/(len-1)
             
            OTOCvalues = OTOCevaluation(W_matrix,V_matrix,time,spectrum)

            times[i] = time
            for j in 1:Dim
                OTOCarrayEnergyScale[i,j] = OTOCvalues[j]

            end
        end

    return OTOCarrayEnergyScale, spectrum, times
end


function GenerateOTOCRegime(N,ξ,ϵ,W,V,nameOTOC;len = 1000,tmin=0,tmax = 10,folder = "OTOCdata/")

    Hamiltonian = Hamiltonian_Nnl(ξ,ϵ,N)

    OTOCarray, spectrum, times = OTOCgeneralRegime(W, V, Hamiltonian,len=len,tmax = tmax,tmin=tmin)

    

    pars = @sprintf " N=%1.0i,chi=%1.2f,eps=%1.2f len = %1.0i" N ξ ϵ len 
    #name = "[l(t),l] abs"
    PathWhole = nameOTOC*pars

    open(folder*PathWhole*"Spectrum.txt", "a") do file
        println(file,spectrum)
    end

    open(folder*PathWhole*"Times.txt", "a") do file
        println(file,times)
    end

    open(folder*PathWhole*"Values.txt", "a") do file
        println(file,OTOCarray)
    end
end

function Read_Results(path::String)

    Data = []
    open(path) do file 
        line = readline(file)
        #line = readline(file)

        line = replace(line, "[" => "")   
        line = replace(line, "]" => "")
        line = replace(line, "," => "")
        line = replace(line, ";" => "")
        line = replace(line, "Any" => "")
        
        elements = split(line)
        Data = [parse(Float64, element) for element in elements]
    end


    return Data
end

function GetOTOCdata(N,ξ,ϵ,len,nameOTOC;folder = "OTOCdata/")

    Dim = Int64((N+2)*(N+1)/2)
    pars = @sprintf " N=%1.0i,chi=%1.2f,eps=%1.2f len = %1.0i" N ξ ϵ len 
    #name = "[l(t),l] abs"
    PathWhole = nameOTOC*pars

    spectrum = Read_Results(folder*PathWhole*"Spectrum.txt")
    OTOCarray = Read_Results(folder*PathWhole*"Values.txt")
    times = Read_Results(folder*PathWhole*"Times.txt")

    #OTOCarray =  reshape(OTOCarray,len,Dim)
    OTOCarray =  reshape(OTOCarray,Dim,len)


    return OTOCarray', spectrum
    
end


    N = 15
    #W2 = W2_Nnl(N)
    #n = n_Nnl(N)
    #l = l_Nnl(N)
   # np = (n+l)/2
  #  nm = (n-l)/2
    Dp = Dp_Nnl(N)
    Dm = Dm_Nnl(N)
 #   Rp = Rp_Nnl(N)
#    Rm = Rm_Nnl(N)

    #O = Dx*Dx

    Dx = (Dp_Nnl(N) + Dm_Nnl(N))/2   
   # Rx = (Rp_Nnl(N) + Rm_Nnl(N))/2   

GenerateOTOCRegime(N,0.1,0.4,Dx,Dx,"[Dx(t),Dx]SHORT",len = 1000,tmin = 0,tmax=100)
GenerateOTOCRegime(N,0.1,0.4,Dx,Dx,"[Dx(t),Dx]LONG",len = 1000,tmin = 1e9,tmax=1e10)
