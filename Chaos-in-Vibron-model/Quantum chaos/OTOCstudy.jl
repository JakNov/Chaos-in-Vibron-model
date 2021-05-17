

#
#
#   zatím nic
#
#





#=
using Plots
using PyCall
using LinearAlgebra


include("SpectrumFunctions.jl")
include("OTOCs.jl")

N = 5



W2 = W2_Nnl(N)
Dp = Dp_Nnl(N)
Dm = Dm_Nnl(N)

n = convert(Array{Float64,2}, n_Nnl(N))
l = convert(Array{Float64,2}, l_Nnl(N))
ns = convert(Array{Float64,2}, ns_Nnl(N))

Rp = Commutator(n,Dp)
Rm = convert(Array{Float64,2},Rp')

W2Cas = 1/2 * (Dp*Dm + Dm*Dp) +l*l

Hamiltonian = n #+ W2
EigSystem = eigen(Hamiltonian)
spectrum = EigSystem.values
vectors = EigSystem.vectors
    
S = vectors
W2_ham = inv(S)*W2*S #n ve vlastni bazi hamiltonianu
n_ham = inv(S)*n*S

W_t(t) = Diagonal(exp.( 1im *t*spectrum)) * W2_ham * Diagonal(exp.( -1im *t*spectrum))
n_t(t) = Diagonal(exp.( 1im *t*spectrum)) * n_ham * Diagonal(exp.( -1im *t*spectrum))

Wt_ex1(t) = Dm*Dp + im*t*(Rp*Dm - Rm*Dp) - sin(t)^2 * (Dm*Dp - Rm*Rp) - im*(t - cos(t)*sin(t))*(Rp*Dm - Rm*Dp) + l + l*l + im*t*(4*ns -2*n)

Wt_ex2(t) = 1/2 *(Dp*Dm +Dm*Dp) + l*l# + im*t/2 * (Rp*Dm + Dm*Rp - Dp*Rm - Rm*Dp) -sin(t)^2 *(Dm*Dp - Rm*Rp) - im*(t - cos(t)*sin(t))*(Rp*Dm - Dp*Rm)



OTOCmatr(t) =  S *(Commutator(n_ham,W_t(t)) * Commutator(n_ham,W_t(t))')* inv(S)

OTOCmatr2(t) =  S *(Commutator(W2_ham,n_t(t)) * Commutator(W2_ham,n_t(t))')* inv(S)


D = Int((N + 1)*(N +2)/2) 

A = Dm*Rp - Rm*Dp
B = Dp*Dm - Rp*Rm

η = Dm*Rp - Rp*Dm

OTOCm2(t) = -cos(2*t)^2 *A*A -im/2 *sin(4*t)*(A*B + B*A) + sin(2*t)^2 * B*B


OTOCm(t) = -Commutator(n,W2)*Commutator(n,W2) + sin(2*t)^2 * (A*A + B*B) - im/2 * sin(4t)*(A*B + B*A) - 2*sin(t)^2 * (A*η + η*A) + im*sin(2t) * (B*η + η*B)

#plotly()



pyplot()
PyPlot.pygui(true)

matr = S*W_t(0.0)*inv(S)
otoc = n*matr - matr*n

otoc2(t) = cos(t)^2 * (n*Dm*Dp - Dm*Dp*n) + sin(t)^2 * (n*Rm*Rp - Rm*Rp*n) + im/2 *sin(2*t) *(n*(Rp*Dm - Rm*Dp) - (Rp*Dm - Rm*Dp)*n)

#a1 = heatmap(convert(Array{Float64,2},real.(matr)))
#a2 = heatmap(convert(Array{Float64,2},imag.(matr)))

a1 = heatmap(convert(Array{Float64,2},real.(Wt_ex1(0.0))))
a2 = heatmap(convert(Array{Float64,2},imag.(Wt_ex1(0.0))))

a1 = heatmap(convert(Array{Float64,2},real.(W2_Nnl(N))))
a2 = heatmap(convert(Array{Float64,2},imag.(W2_Nnl(N))))

a3 = heatmap(convert(Array{Float64,2},real.(Wt_ex2(0.0))))
a4 = heatmap(convert(Array{Float64,2},imag.(Wt_ex2(0.0))))

a1 = heatmap(convert(Array{Float64,2},real.(W2*Dp -Dp*W2)))
a2 = heatmap(convert(Array{Float64,2},imag.(W2*Dp -Dp*W2)))

a3 = heatmap(convert(Array{Float64,2},real.(W2Cas*Dp -Dp*W2Cas)))
a4 = heatmap(convert(Array{Float64,2},imag.(W2Cas*Dp -Dp*W2Cas)))

a1 = heatmap(convert(Array{Float64,2},real.(W2)))
a2 = heatmap(convert(Array{Float64,2},imag.(W2)))

a1 = heatmap(convert(Array{Float64,2},real.(W2Cas)))
a2 = heatmap(convert(Array{Float64,2},imag.(W2Cas)))

a1 = heatmap(convert(Array{Float64,2},real.(Dx_Nnl(N))))
a2 = heatmap(convert(Array{Float64,2},imag.(Dx_Nnl(N))))

#a3 = heatmap(convert(Array{Float64,2},real.(1/2 * (Dp*Dm) )))
#a4 = heatmap(convert(Array{Float64,2},imag.(1/2 * (Dp*Dm) )))

a3 = heatmap(convert(Array{Float64,2},real.(1/2*(Dp + Dm))))
a4 = heatmap(convert(Array{Float64,2},imag.(1/2*(Dp + Dm))))

a5 = heatmap(convert(Array{Float64,2},-2*real.(Hamiltonian_Nnl(0.0,1.0,N) - n)))
a6 = heatmap(convert(Array{Float64,2},-2*imag.(Hamiltonian_Nnl(0.0,1.0,N) - n)))
plot(a1,a3,a5)

WW = 1/2 * (Dp*Dm + Dm*Dp) + l*l

a1 = heatmap(convert(Array{Float64,2},real.(Dp*Dm - Dm*Dp)))
a2 = heatmap(convert(Array{Float64,2},imag.(Dp*Dm - Dm*Dp)))

#a3 = heatmap(convert(Array{Float64,2},real.(Hamiltonian_Nnl(1.0,0.0,N))))
#a4 = heatmap(convert(Array{Float64,2},imag.(Hamiltonian_Nnl(1.0,0.0,N))))



a3 = heatmap(convert(Array{Float64,2},real.(W2*Dp - Dp*W2)))
a4 = heatmap(convert(Array{Float64,2},imag.(W2*Dp - Dp*W2)))

a1 = heatmap(convert(Array{Float64,2},real.(WW*Dp - Dp*WW - (WW*Dm - Dm*WW))))
a2 = heatmap(convert(Array{Float64,2},imag.(WW*Dp - Dp*WW- (WW*Dm - Dm*WW))))

a3 = heatmap(convert(Array{Float64,2},real.(WW*Dm - Dm*WW)))
a4 = heatmap(convert(Array{Float64,2},imag.(WW*Dm - Dm*WW)))

Dm2 = Dm_Nnl(N)

WWW = 1/2 * (Dp*Dm2 + Dm2*Dp) + l*l

a1 = heatmap(convert(Array{Float64,2},real.(Dp)))
a2 = heatmap(convert(Array{Float64,2},imag.(Dp)))

a3 = heatmap(convert(Array{Float64,2},real.(Dm2)))
a4 = heatmap(convert(Array{Float64,2},imag.(Dm2)))

plot(a1,a2,a3,a4)
#=
matrix1 = real.(OTOCmatr(0.0))
matrix2 = real.(OTOCmatr2(0.0))
#p1 = heatmap(convert(Array{Float64,2},real.(OTOCm(10.0)*OTOCm(10.0)')))
p2 = heatmap(convert(Array{Float64,2},matrix1))
p3 = heatmap(convert(Array{Float64,2},matrix2))
println(tr(matrix1))
println(tr(matrix2))

#p4 = heatmap(convert(Array{Float64,2},imag.(OTOCmatr(1.0))))
#p5 = heatmap(convert(Array{Float64,2},imag.(OTOCmatr2(1.0))))

plot(p2,p3)
#heatmap(convert(Array{Float64,2},Commutator(n,W2)))
plot(eigen(matrix1).values)
plot!(eigen(matrix2).values)
plot!(eigen(real.(OTOCmatr(1.0))).values)
plot!(eigen(real.(OTOCmatr2(1.0))).values)
=#
=#