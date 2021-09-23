using Plots
using PyCall
using Printf
using LinearAlgebra

#
#
#   Jak vypadají vlastní báze hamiltoniánů - míchání podprostorů n a j 
#
#

include("SpectrumFunctions.jl")
include("ChaosIndicators.jl")


ξ = 0.5
ϵ = 0.2

N = 30
D =  Int((N + 1)*(N +2)/2) 

Hamiltonian = Hamiltonian_Nnl(ξ,ϵ,N)
EigSystem = eigen(Hamiltonian)

spectrum = EigSystem.values
vectors = EigSystem.vectors
S2 = vectors


plotly()
for i in 1:length(spectrum)
    E = spectrum[i]
    p1 = scatter(vectors[:,i].*vectors[:,i], label = "level num $i")
    p1 = plot!(title = "E = $E")

    Name_fig = "vector_basis chi = $ξ eps = $ϵ N = $N level $i"
    plot(p1)
    plot!(size=(500,500))
    savefig(Name_fig )
end






