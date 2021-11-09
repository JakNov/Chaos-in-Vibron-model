using Plots
using PyCall
using Printf
using LinearAlgebra


include("SpectrumFunctions.jl")
include("ChaosIndicators.jl")
include("OTOCs.jl")

#
#
#       Vibronový a BEC hamiltonián  - numericky získaná matice přechodu - problémy
#
#

N = 30
ξ = 0.2
ϵ = 0.0

D =Int((N+1)*(N+2)/2)
Num_0 = Int((N - N%2)/2 + 1)

S = TrMaNlnTolml(N)

Hr = HamiltonianR_Nln(ξ,ϵ,N)
Hd = Hamiltonian_Nln(ξ,ϵ,N)



Ud = eigen(Hd).vectors
U = eigen(inv(S)*Hr*S).vectors
Udd = Ud*inv(U)

BlokDiagonalHamiltonianD = inv(Udd)*Hd*Udd

BlokDiagonalHamiltonianD1 = BlokDiagonalHamiltonianD[1:Int((D - Num_0)/2) + Num_0,1:Int((D - Num_0)/2) + Num_0]
BlokDiagonalHamiltonianD2 = BlokDiagonalHamiltonianD[Int((D - Num_0)/2) + Num_0 + 1:end,Int((D - Num_0)/2) + Num_0 + 1:end]

system1 = eigen(BlokDiagonalHamiltonianD1)
system2 = eigen(BlokDiagonalHamiltonianD2)

#=
#
#   Blok na spektrlní statistiky - nemá smysl pro malá N
#

nnd1 = NND(real.(system1.values))
nnd2 = NND(real.(system2.values))

b1 = BrodyParameter(nnd1)
b2 = BrodyParameter(nnd2)



pyplot()
PyPlot.pygui(true)

#p1 = heatmap(BlokDiagonalHamiltonianD)
#p2 = heatmap(BlokDiagonalHamiltonianD1)
#p3 = heatmap(BlokDiagonalHamiltonianD2)

d1 = plot(nnd1)
d1 = plot!([i/10 for i in 1:100],[Poisson(i/10) for i in 1:100], label = "Poisson")
d1 = plot!([i/10 for i in 1:100],[Brody(i/10, b1) for i in 1:100], label = "Brody")

d2 = plot(nnd2)
d2 = plot!([i/10 for i in 1:100],[Poisson(i/10) for i in 1:100], label = "Poisson")
d2 = plot!([i/10 for i in 1:100],[Brody(i/10, b2) for i in 1:100], label = "Brody")

plot(d1,d2)

println(b1)
println(b2)

=#

    #
    # OTOCy v bázi, ve které je hamiltonián blokově diagonální
    #
    
    TrM = TrMaNnlToNln(N)

    W = inv(Udd)*inv(TrM)*W2_Nnl(N)*TrM*Udd #blokove diagonalni baze
    V = inv(Udd)*n_Nln(N)*Udd



    W1 = W[1:Int((D - Num_0)/2) + Num_0,1:Int((D - Num_0)/2) + Num_0]
    W2 = W[Int((D - Num_0)/2) + Num_0 + 1:end,Int((D - Num_0)/2) + Num_0 + 1:end]


    V1 = V[1:Int((D - Num_0)/2) + Num_0,1:Int((D - Num_0)/2) + Num_0]
    V2 = V[Int((D - Num_0)/2) + Num_0 + 1:end,Int((D - Num_0)/2) + Num_0 + 1:end]
    
    S1 = system1.vectors
    S2 = system2.vectors

    W1_matrix = inv(S1)*W1*S1 #n ve vlastni bazi hamiltonianu
    V1_matrix = inv(S1)*V1*S1
    
    W2_matrix = inv(S2)*W2*S2 #n ve vlastni bazi hamiltonianu
    V2_matrix = inv(S2)*V2*S2

    M = zeros(Float64,D,D)
    M[1:Int((D - Num_0)/2) + Num_0,1:Int((D - Num_0)/2) + Num_0] = S1
    M[Int((D - Num_0)/2) + Num_0 + 1:end,Int((D - Num_0)/2) + Num_0 + 1:end] = S2


    system = eigen(BlokDiagonalHamiltonianD)
    #Ss = system.vectors

    W_matrix = inv(M)*W*M
    V_matrix = inv(M)*V*M
    
      
    len = 1000

    OTOCarray1 = Array{Float64,1}(undef, len)
    OTOCarray2 = Array{Float64,1}(undef, len)
    OTOCarray0 = Array{Float64,1}(undef, len)



    for i in 1:len
        println(i)
        OTOCarray1[i] = real(OTOC(W1_matrix,V1_matrix,2*(i-1)/20,system1.values)*length(system1.values)/D)
        OTOCarray2[i] = real(OTOC(W2_matrix,V2_matrix,2*(i-1)/20,system2.values)*length(system2.values)/D)
        OTOCarray0[i] = real(OTOC(W_matrix,V_matrix,2*(i-1)/20,vcat(system1.values,system2.values)))

    end

    OTOCarrayWhole = OTOCarray(ξ,ϵ,N,len=len, Vmatrix = "W2", Wmatrix = "n")
    #OTOCarrayWhole2 = OTOCarray(ξ,ϵ,N,len=len, Vmatrix = "n", Wmatrix = "W2")

    pyplot()
    PyPlot.pygui(true)

    plot(OTOCarray1)
    plot!(OTOCarray2)
    plot!(OTOCarray0)

    plot!(OTOCarrayWhole)
    plot!(OTOCarray1 + OTOCarray2)
   # plot!(OTOCarrayWhole2)



