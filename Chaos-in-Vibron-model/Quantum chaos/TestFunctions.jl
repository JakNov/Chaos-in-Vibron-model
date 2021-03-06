using Plots
using PyCall
using Printf
using SparseArrays
using ArnoldiMethod
using Arpack

include("SpectrumFunctions.jl")
include("ChaosIndicators.jl")

#
#   Plotting - ξ dependence of levels
#

Val_ξ = Float64[]
spc=Float64[]
N = 20
ϵ = 0.0

n = 75
for i in 1:n
    ξ = (i-1)/(n-1)
    append!(Val_ξ,ξ)
    spectrum = eigen(Hamiltonian_Nln(ξ,ϵ,N)).values
    append!(spc,spectrum)
end

spc = transpose(reshape(spc,Int((N + 1)*(N +2)/2),n))

pyplot()
PyPlot.pygui(true)
plot(Val_ξ,spc, title = "Přechod v ξ",legend = false, c = "black")
xlabel!("ξ")
ylabel!("Energie")


#
#   Hamiltonian plot heat map
#



pyplot()
PyPlot.pygui(true)
heatmap(Hamiltonian_Nln(0.5,0.5,10),c = :balance, clim = (-10,10))




#
#  
#
#=
B_array = Float64[]
σ_array = Float64[]
x_array = Float64[]
η_array = Float64[]
ξ_E_array = Float64[]

N = 50
for i in 1:2N

    println(i)
    #Hamiltonian = Hamiltonian_Nln((i-1)/N,0.5,75)
    Hamiltonian = Hamiltonian_Nln(0.5,(i-1)/N,75)
    EigSystem = eigen(Hamiltonian)

    spectrum = EigSystem.values
    vectors = EigSystem.vectors
    #spectrum = CutDegeneracy(EigSystem.values) # co s ni?

    HistSpacing = NND(spectrum)
    #plot(HistSpacing)
    B,σ  = BrodyParameter(HistSpacing)
    Eta = η(spectrum)
    ξ_E = IPR(vectors)

    plot(HistSpacing, label = "NND")
    xdata = LinRange(0.0,5.0,100)
    P = Poisson.(xdata)
    plot!(xdata,P,label = "Poisson")
    W = Wigner.(xdata)
    plot!(xdata,W,label = "Wigner")
    Br_y = Brody(xdata,B)
    name = @sprintf "Brody %1.3f +/- %1.3f" B σ
    plot!(xdata,Br_y,label = name)
    Name_fig = @sprintf "Test_eps=%1.2fxi=%1.2fN=%3.0i" 0.5 (i-1)/N 75
    savefig(Name_fig)

    append!(B_array,B)
    append!(σ_array,σ)
    append!(η_array,Eta)
    append!(ξ_E_array,ξ_E)
    append!(x_array,(i-1)/N)
end

pyplot()
plot(x_array,B_array,yerror=σ_array,title = "ξ = 0.5",label = "Brody",xlabel = "ϵ", ylabel = "β")
plot!(x_array,η_array,label = "η")
plot!(x_array,ξ_E_array,label = "ξ_E")
#savefig("")
=#