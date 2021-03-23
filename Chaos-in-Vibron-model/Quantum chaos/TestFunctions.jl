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
#=
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

=#


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

#
#   Plotting - vectors and spectrum
#

#=
ξ = 0.05
ϵ = 0.05

N =  8#10

Hamiltonian = Hamiltonian_Nln(ξ,ϵ,N)
EigSystem = eigen(Hamiltonian)

spectrum = EigSystem.values
vectors = EigSystem.vectors

println(length(spectrum))
println(length(CutDegeneracy(spectrum)))

plotly()
#PyPlot.pygui(true)
heatmap(vectors, c= :blues)
scatter!(spectrum)
scatter!(CutDegeneracy(spectrum),label = "cut")

num = 0
if N%2 == 0
    
    for j in 1:Int64(N/2)
        plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + j for i in 1:Int64((N+1)*(N+2)/2)],c = :black)
        global num = num + j
        plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + j for i in 1:Int64((N+1)*(N+2)/2)],c = :black)
        global num = num + j
    end
    
    plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + N/2+1 for i in 1:Int64((N+1)*(N+2)/2)],c = :black)

end 
Name_fig = "vectors0.05 0.05 8"
savefig(Name_fig)
=#



#
#
#

#=
ξ = 0.5
ϵ = 0.0

Hamiltonian = Hamiltonian_Nnl(ξ,ϵ,20)
EigSystem = eigen(Hamiltonian)

spectrum = EigSystem.values
vectors = EigSystem.vectors


pyplot()
PyPlot.pygui(true)
heatmap(vectors, c= :blues)
scatter!(spectrum)=#
#=
ξ = 0.05
ϵ = 0.05

N =  8#10

Hamiltonian = Hamiltonian_Nln(ξ,ϵ,N)
EigSystem = eigen(Hamiltonian)

spectrum = EigSystem.values
vectors = EigSystem.vectors
v1 = vectors[:,37]
v2 = vectors[:,38]
len = length(v2)
v3 = [v2[len - i+1] for i in 1:len]


pyplot()
PyPlot.pygui(true)
scatter(v1,c=:red)
scatter!(v2,c=:blue)
plot!(v1,c=:red)
plot!(v2,c=:blue)
=#
#=
ξ = 0.05
ϵ = 0.0

N =  8#10

Hamiltonian = Hamiltonian_Nln(ξ,ϵ,N)
EigSystem = eigen(Hamiltonian)

spectrum = EigSystem.values
vectors = EigSystem.vectors
v1 = vectors[:,37]
v2 = vectors[:,38]
len = length(v2)
v3 = [v2[len - i+1] for i in 1:len]


pyplot()
PyPlot.pygui(true)
scatter!(v1)
scatter!(v3)

=#
#
#   Plotting - Mean values - l, n
#
#

#=
ξ = 0.3
ϵ = 0.0#01

N =  75

Hamiltonian = Hamiltonian_Nln(ξ,ϵ,N)
EigSystem = eigen(Hamiltonian)

spectrum = EigSystem.values
vectors = EigSystem.vectors

len =length(spectrum)
amplitudes = vectors .* vectors
l = []
n = []

Basis = Basis_Nln(N)

for i in 1:len
    append!(l,dot(amplitudes[:,i],Basis[2,:]))
    append!(n,dot(amplitudes[:,i],Basis[3,:]))
end
#plotly()
pyplot()
PyPlot.pygui(true)
title = @sprintf "Střední hodnoty vlastních stavů \n N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
p1 = scatter(spectrum,n, xlabel= "Energie", ylabel = "n", title = title, label = "<n>")
p2 = scatter(spectrum,l, xlabel= "Energie", ylabel = "l", label = "<l>")

plot(p1, p2, layout = (1,2))
=#

#
#   NNspacing - l,n subspaces
#

#=
N = 75
ξ = 0.5
ϵ = 0.0

Hamiltonian = Hamiltonian_Nln(ξ,ϵ,N)
EigSystem = eigen(Hamiltonian)
spectrum = EigSystem.values
vectors = EigSystem.vectors

Basis = Basis_Nln(N)
len = length(spectrum)

SpectrumPositive_ln1 = Float64[]
SpectrumPositive_ln2 = Float64[]

SpectrumNegative_ln1 = Float64[]
SpectrumNegative_ln2 = Float64[]

for i in 1:len
    
    if Basis[2,i] >= 0 && Basis[3,i]%2 != 0
        append!(SpectrumPositive_ln1, spectrum[i])
    end
    if Basis[2,i] >= 0 && Basis[3,i]%2 == 0
        append!(SpectrumPositive_ln2, spectrum[i])
    end

    if Basis[2,i] <= 0 && Basis[3,i]%2 != 0
        append!(SpectrumNegative_ln1, spectrum[i])
    end

    if Basis[2,i] <= 0 && Basis[3,i]%2 == 0
        append!(SpectrumNegative_ln2, spectrum[i])
    end
end


HistSpacing_ln1 = NND(SpectrumPositive_ln1)
HistSpacing_ln2 = NND(SpectrumPositive_ln2)

HistSpacing_nln1 = NND(SpectrumNegative_ln1)
HistSpacing_nln2 = NND(SpectrumNegative_ln2)


pyplot()
#PyPlot.pygui(true)

xdata = LinRange(0.0,5.0,100)
P = Poisson.(xdata)

p1 = plot(HistSpacing_ln1, label = "NNS",title = "Kladne l, n liche")
p1 = plot!(xdata,P,label = "Poisson")
p2 = plot(HistSpacing_ln2, label = "NNS",title = "Kladne l, n sude")
p2 = plot!(xdata,P,label = "Poisson")
p3 = plot(HistSpacing_nln1, label = "NNS",title = "Zaporne l, n liche")
p3 = plot!(xdata,P,label = "Poisson")
p4 = plot(HistSpacing_nln2, label = "NNS",title = "Zaporne l, n sude")
p4 = plot!(xdata,P,label = "Poisson")
p5 = plot(NND(spectrum), label = "NNS", title = "Spektrum")
p5 = plot!(xdata,P,label = "Poisson")
p6 = plot(spectrum)




plot(p1, p2, p3, p4 ,p5 ,p6, layout = (2, 3), legend = false)
=#