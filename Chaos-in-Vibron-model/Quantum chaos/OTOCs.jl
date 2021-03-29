using Plots
using PyCall
using Printf
using SparseArrays
using LinearAlgebra
using StatsBase
using FFTW

include("SpectrumFunctions.jl")

#
#   OTOC  C(t) = <[W(t),V]^deg[W(t),V]>
#
function OTOC(W::Array{Float64,2},V::Array{Float64,2},t::Float64,spectrum::Array{Float64,1})

    Dim = length(spectrum)
    W_t = Diagonal(exp.( 1im *t*spectrum)) * W * Diagonal(exp.( -1im *t*spectrum))

    WV = W_t * V
    VW = V * W_t

    OTOC = 2/Dim * (tr((VW - WV)* WV))
    return OTOC
end

#
#   OTOC - evaluation of OTOC C(0 -> len, step = len/50) = <[W(t),V]^deg[W(t),V]>
#

function OTOCarray(ξ::Float64,ϵ::Float64,N::Int64; len::Int64 = 10000, basis::String = "Nnl", Vmatrix = "W2", Wmatrix = "n")

    if basis == "Nnl"
    
        Hamiltonian = Hamiltonian_Nnl(ξ,ϵ,N)
            
        if Vmatrix == "W2"
            V = W2_Nnl(N)
        elseif Vmatrix == "n"
            V = n_Nnl(N)
        elseif Vmatrix == "Dx"
            V = Dx_Nnl(N)
        end

        if Wmatrix == "W2"
            W = W2_Nnl(N)
        elseif Wmatrix == "n"
            W = n_Nnl(N)
        elseif Wmatrix == "Dx"
            W = Dx_Nnl(N)
        end
            
    elseif basis == "Nln"
    
        Hamiltonian = Hamiltonian_Nln(ξ,ϵ,N)
    
        if Vmatrix == "W2"
            V = W2_Nln(N)
        elseif Vmatrix == "n"
            V = n_Nln(N)
        end
    
        if Wmatrix == "W2"
             W = W2_Nln(N)
        elseif Wmatrix == "n"
            W = n_Nln(N)
        end
    
    else 
        println("UFF")
    end
    
    
        EigSystem = eigen(Hamiltonian)
        spectrum = EigSystem.values
        vectors = EigSystem.vectors
    
        Dim = length(spectrum)
    
        S = vectors
        W_matrix = inv(S)*W*S #n ve vlastni bazi hamiltonianu
        V_matrix = inv(S)*V*S
    
    
            #@time OTOCarray = [real(OTOC(n_matrix,W2_matrix,1.0*((i-1)/(50-1)),spectrum)) for i in 1:len]
    
        OTOCarray = Array{Float64,1}(undef, len)

        Threads.@threads for i in 1:len
            OTOCarray[i] = real(OTOC(W_matrix,V_matrix,1.0*((i-1)/(50-1)),spectrum))

            #println("Point $i on thread $(Threads.threadid())")

        end
    
    
    return OTOCarray
    
    
end




#=
n_t2 = ComplexF64[]


@time for i in 1:Dim
    for j in 1:Dim
        append!(n_t2,n_matrix[i,j]*exp(1im*t*(spectrum[i]-spectrum[j])))
    end
end

n_t2 = reshape(n_t2,(Dim,:));=#
