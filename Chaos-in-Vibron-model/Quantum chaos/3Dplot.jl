using Plots
using PyCall

include("SpectrumFunctions.jl")

#
#       plotting the wavefunction
#

N = 20
dim = Int64((N+1)*(N+2)/2)
Qp = Qp_Nnl(N)
Qm = Qm_Nnl(N)

Qa = -1/sqrt(2) *(Qp - Qm)
Qb = im/sqrt(2) *(Qp + Qm)

SystemA = eigen(Qa)
SystemB = eigen(Qb)

ξ = 0.2
ϵ = 0.0
H = Hamiltonian_Nnl(ξ,ϵ,N)
SystemH = eigen(H)

i = 1 #hladina
ψa = [real(SystemA.vectors[:,j]'*SystemH.vectors[:,i]) for j in 1:dim]
ψb = [real(SystemB.vectors[:,j]'*SystemH.vectors[:,i]) for j in 1:dim]

Sa = ψa.* conj.(ψa)
Sb = ψb.* conj.(ψb)

WaveFunction = Sa .* real.(Sb)'

n = 35#Int64((dim + dim%2)/2)

x = SystemA.values
y = SystemB.values

xmin = x[1]
xmax = x[end]
xstep = (xmax-xmin)/n

ymin = y[1]
ymax = y[end]
ystep = (ymax-ymin)/n

indx = 1
IndX = [1]

indy = 1
IndY = [1]

for j in 2:dim
   # println([x[j]],[xmin + xstep*(indx)],[x[j-1]])
    if x[j] > (xmin + xstep*(indx)) > x[j-1]
        append!(IndX,[j-1])
        global indx=indx+1
        #println(indx)
    end

    if y[j] > (ymin + ystep*(indy)) > y[j-1]
        append!(IndY,[j-1])
        global indy=indy+1
    end

end

append!(IndX,[dim])
append!(IndY,[dim])


lenx = length(IndX)
leny = length(IndY)

Z = zeros(Float64,lenx,leny)

for j in 1:lenx-1
    for k in 1:leny-1
        Z[j,k] = sum(WaveFunction[IndX[j]:IndX[j+1],IndY[k]:IndY[k+1]])
    end
end

#=
x = []
y = []
z = []

for j in 1:dim

    for k in 1:dim
        append!(x,[ψa[j]])
        append!(y,[real(ψb[k])])
        append!(z,[WaveFunction[j,k]])    
    end
end
=#

pyplot()
PyPlot.pygui(true)
heatmap(Z)
#surface(SystemA.values,SystemB.values,Sa .* real.(Sb)')
#scatter(x,y,z)