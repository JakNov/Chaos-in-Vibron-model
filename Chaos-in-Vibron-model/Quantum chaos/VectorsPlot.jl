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


ξ = 0.05
ϵ = 0.05

N = 4#10
D =  Int((N + 1)*(N +2)/2) 

HamiltonianR = HamiltonianR_Nln(ξ,ϵ,N)
EigSystemR = eigen(HamiltonianR)

spectrumR = EigSystemR.values
vectorsR = EigSystemR.vectors
S2 = vectorsR

HamiltonianD = Hamiltonian_Nln(ξ,ϵ,N)
EigSystemD = eigen(HamiltonianD)

spectrumD = EigSystemD.values
vectorsD = EigSystemD.vectors
S1 = vectorsD

#
#   operatory
#

n = n_Nnl(N)
l = l_Nnl(N)

Dp = Dp_Nnl(N)
Rp = Rp_Nnl(N)

HD = #=(1 - ξ)*n - ξ/(N-1)*(1/2 * (Dp*Dp' + Dp'*Dp) + l*l) + ϵ/2 *=# (Dp + Dp')
HR = #=(1 - ξ)*n - ξ/(N-1)*(1/2 * (Rp*Rp' + Rp'*Rp)+ l*l) + ϵ/2 *=# (Rp + Rp')

HDL =S2* inv(S1) * HD * S1*inv(S2)
HRL = S2*inv(S1) * HR * S1*inv(S2)


plotly()
#PyPlot.pygui(true)
#p1 = heatmap(eigen(HDL).vectors, c= :blues)

#vectors = eigen(HRL).vectors
#p1 = heatmap(vectorsD-vectorsR, c= :blues)
p1 = heatmap(vectorsD, c= :blues)

#p1 = scatter!(spectrumL)

num = 0
if N%2 == 0
    
    for j in 1:Int64(N/2)
        p1 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + j + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :black,label = false)
        global num = num + j
        p1 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + j + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :black,label = false)
        global num = num + j
    end
    
    #plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + N/2+1 for i in 1:Int64((N+1)*(N+2)/2)],c = :black)

end 

global num = 0
if N%2 == 0
    
    for j in 1:Int64(N/2)
        p1 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[D- (num + j) + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :red,label = false)
        global num = num + j
        p1 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[D - (num + j) + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :red,label = false)
        global num = num + j
    end
end 

p1 = plot!(title = "D hamiltonian")

#p2 = heatmap(vectorsD, c= :blues)
#p2 = heatmap(Hamiltonian_Nln(ξ,ϵ,N), c= :blues)
#p2 = heatmap(vectorsD + vectorsR, c= :blues)
p2 = heatmap(vectorsR, c= :blues)


#scatter!(spectrum)

num = 0
if N%2 == 0
    
    for j in 1:Int64(N/2)
        p2 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + j + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :black,label = false)
        global num = num + j
        p2 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + j + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :black,label = false)
        global num = num + j
    end
    
    #plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[num + N/2+1 for i in 1:Int64((N+1)*(N+2)/2)],c = :black)

end 

global num = 0
if N%2 == 0
    
    for j in 1:Int64(N/2)
        p2 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[D- (num + j) + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :red,label = false)
        global num = num + j
        p2 = plot!([i for i in 1:Int64((N+1)*(N+2)/2)],[D - (num + j) + 0.5 for i in 1:Int64((N+1)*(N+2)/2)],c = :red,label = false)
        global num = num + j
    end
end 

p2 = plot!(title = "R hamiltonian")


Name_fig = "vector_basis"
plot(p1,p2)
plot!(size=(2000,1000))
savefig(Name_fig)

