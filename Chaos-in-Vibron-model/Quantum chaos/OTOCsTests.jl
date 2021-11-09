using Plots
using PyCall
using Printf
using StatsBase
using LinearAlgebra


include("SpectrumFunctions.jl")
include("OTOCs.jl")

#
#
#   Plot více OTOCů zároveň
#
#

ξ = 0.001
ϵ = 0.0
N = 15

#
# pred spustenim julie 1.3 export do terminalu ~$ JULIA_NUM_THREADS=15
#

println("Num of threads: $(Threads.nthreads())") #number of threads 


len = Int64(1e3) #time


@time otoc_nn = OTOCarray(ξ,ϵ,N,Vmatrix = "n",Wmatrix = "n", len = len)
@time otoc_Wn = OTOCarray(ξ,ϵ,N,Vmatrix = "W2",Wmatrix = "n", len = len)
@time otoc_WW = OTOCarray(ξ,ϵ,N,Vmatrix = "W2",Wmatrix = "W2", len = len)
@time otoc_nW = OTOCarray(ξ,ϵ,N,Vmatrix = "n",Wmatrix = "W2", len = len)
@time otoc_DxW = OTOCarray(ξ,ϵ,N,Vmatrix = "Dx",Wmatrix = "W2", len = len)
@time otoc_WDx = OTOCarray(ξ,ϵ,N,Vmatrix = "W2",Wmatrix = "Dx", len = len)
@time otoc_DxDx = OTOCarray(ξ,ϵ,N,Vmatrix = "Dx",Wmatrix = "Dx", len = len)
@time otoc_nDx = OTOCarray(ξ,ϵ,N,Vmatrix = "n",Wmatrix = "Dx", len = len)
@time otoc_Dxn = OTOCarray(ξ,ϵ,N,Vmatrix = "Dx",Wmatrix = "n", len = len)


title = @sprintf "OTOC"
label1 =  @sprintf "OTOC [n,n(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label2 =  @sprintf "OTOC [W^2,n(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label3 =  @sprintf "OTOC [W^2,W^2(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label4 =  @sprintf "OTOC [n,W^2(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label5 =  @sprintf "OTOC [Dx,W^2(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label6 =  @sprintf "OTOC [Dx,Dx(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label7 =  @sprintf "OTOC [n,Dx(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label8 =  @sprintf "OTOC [W^2,Dx(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label9 =  @sprintf "OTOC [Dx,n(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ



plotly()


otoc_nn = otoc_nn * 1/(mean(otoc_nn) + 1)
otoc_Wn = otoc_Wn * 1/(mean(otoc_Wn) + 1)
otoc_WW = otoc_WW * 1/(mean(otoc_WW) + 1)
otoc_nW = otoc_nW * 1/(mean(otoc_nW) + 1)
otoc_DxW = otoc_DxW * 1/(mean(otoc_DxW) + 1)
otoc_WDx = otoc_WDx * 1/(mean(otoc_WDx) + 1)
otoc_DxDx = otoc_DxDx * 1/(mean(otoc_DxDx) + 1)
otoc_nDx = otoc_nDx * 1/(mean(otoc_nDx) + 1)
otoc_Dxn = otoc_Dxn * 1/(mean(otoc_Dxn) + 1)



plot(otoc_nn, title = title,label = label1, size = (1500, 800),xlabel = "time",ylabel = "C(t)/<C(t) + 1>")
plot!(otoc_Wn, label = label2)
plot!(otoc_WW, label = label3)
plot!(otoc_nW, label = label4)
plot!(otoc_DxW, label = label5)
plot!(otoc_WDx, label = label8)
plot!(otoc_DxDx, label = label6)
plot!(otoc_nDx, label = label7)
plot!(otoc_Dxn, label = label9)

nameplot =  @sprintf "OTOC N=%3.0i ξ = %1.2f ϵ = %1.3f length = 5e4" N ξ ϵ

savefig(nameplot)



#
#   Fast Fourier Transform
#
#= 
otoc_plot = plot(otoc_nn,label = label1)

F = fft(otoc_nn) |> fftshift
freqs = fftfreq(len, 1.0*49.0) |> fftshift

otoc_freqs = plot(freqs,abs.(F),label = "Frekvence")
plot(otoc_plot, otoc_freqs, layout = 2)
=#


#
#   Zatim nefunkcni
#


#=
ξ = 0.0
ϵ = 0.0
N = 5

#
# pred spustenim julie 1.3 export do terminalu ~$ JULIA_NUM_THREADS=15
#

println("Num of threads: $(Threads.nthreads())") #number of threads 


len = Int64(5000) #time


@time otoc_Wn = OTOCarray(ξ,ϵ,N,Vmatrix = "W2",Wmatrix = "n", len = len)
@time otoc_nW = OTOCarray(ξ,ϵ,N,Vmatrix = "n",Wmatrix = "W2", len = len)
title = @sprintf "OTOC"
label =  @sprintf "OTOC [W^2,n(t)]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ
label2 =  @sprintf "OTOC [n(t),W^2]"# N=%3.0i ξ = %1.2f ϵ = %1.3f " N ξ ϵ


W2 = W2_Nnl(N)
Dp = Dp_Nnl(N)
Dm = convert(Array{Float64,2},Dp')

n = convert(Array{Float64,2}, n_Nnl(N))
l = convert(Array{Float64,2}, l_Nnl(N))
ns = convert(Array{Float64,2}, ns_Nnl(N))

Rp = -Commutator(n,Dp)
Rm = convert(Array{Float64,2},Rp')

f1(t) = - cos(t)*cosh(t)
f2(t) = - im*(cosh(t)*sin(t) + cos(t)*sinh(t))
f3(t) = - sin(t)*sinh(t)
f4(t) = im*(cosh(t)*sin(t) - cos(t)*sinh(t))

OT(t) = f1(t)*(Rp*Dm + Rm*Dp) + f2(t)*(Rp*Rm + l) + f3(t)*(Rm*Dp -Rp*Dm + 4*ns -2*n) + f4(t)*(Dp*Dm + l)
TO(t) = f1(t)*(Dm*Rp + Dp*Rm) - f2(t)*(Rp*Rm + l) + f3(t)*(Dm*Rp -Dp*Rm + 4*ns -2*n) - f4(t)*(Dp*Dm + l)

OTOCapproximated = Array{Float64,1}(undef, len)

D = Int((N + 1)*(N +2)/2) 

A = Dm*Rp - Rm*Dp
B = Dp*Dm - Rp*Rm

η = Dm*Rp - Rp*Dm

C(t) = -tr(Commutator(n,W2)*Commutator(n,W2)) + sin(2*t)^2 * tr(A*A + B*B) - im/2 * sin(4t)*tr(A*B + B*A) - 2*sin(t)^2 * tr(A*η + η*A) + im*sin(2t) * tr(B*η + η*B)
Coperator(t) = #=-Commutator(n,W2)*Commutator(n,W2) +=# sin(2*t)^2 * (A*A + B*B) - im/2 * sin(4t)*(A*B + B*A) - 2*sin(t)^2 * (A*η + η*A) + im*sin(2t) * (B*η + η*B)

Creal(t) = -D^2 + sin(2*t)^2 *D^2 - 1/2 * sin(4t)*D^2 - 2*sin(t)^2 *D^2 + 1*sin(2t) *D^2




#=function Cfunction(t)
   
    D = Int((N + 1)*(N +2)/2)
    
    = -tr(Commutator(n,W2)*Commutator(n,W2)) + sin(2*t)^2 * tr(A*A + B*B) - im/2 * sin(4t)*tr(A*B + B*A) - 2*sin(t)^2 * tr(A*η + η*A) + im*sin(2t) * tr(B*η + η*B)

end=#

for i in 1:len
    OTOCapproximated[i] = real(C(i)/D)
    println(real(C(i)/D))
    println(imag(C(i)/D))
    #println(i)
end


#=
f1(t) = - im*(cosh(t/sqrt(2))*sin(t/sqrt(2)) + cos(t/sqrt(2))*sinh(t/sqrt(2)))/(sqrt(2))
f2(t) = sin(t/sqrt(2))*sinh(t/sqrt(2))
f3(t) = im*(cosh(t/sqrt(2))*sin(t/sqrt(2)) - cos(t/sqrt(2))*sinh(t/sqrt(2)))/(sqrt(2))
f4(t) = -1 + cos(t/sqrt(2))*cosh(t/sqrt(2))


W2 = W2_Nnl(N)
Dp = Dp_Nnl(N)
Dm = convert(Array{Float64,2},Dp')

n = convert(Array{Float64,2}, n_Nnl(N))

Rp = -Commutator(n,Dp)
Rm = convert(Array{Float64,2},Rp')

A1 = 1/2 * (AntiCommutator(Rp,Dm) + AntiCommutator(Rm,Dp))
A2 = -AntiCommutator(Rp,Rm)
A3 = -(-AntiCommutator(Rp,Dm) + AntiCommutator(Rm,Dp))
A4 = 2*AntiCommutator(Dp,Dm)

nt(t) = n + ξ*(f1(t)*A1 + f2(t)*A2 + f3(t)*A3 + f4(t)*A4)


function OTOCapprox(W,V)
    Dim = Int((N + 1)*(N +2)/2)
    OTOC = 2/Dim * (tr(Commutator(W,V) * (Commutator(W,V)')))
    return OTOC
end

OTOCapproximated = Array{Float64,1}(undef, len)

for i in 1:len
    OTOCapproximated[i] = real(OTOCapprox(nt(((i-1)/(50-1))),W2))
    println(i)
end
=#


plotly()
plot(otoc_Wn, title = title,label = label, size = (1500, 800),xlabel = "time",ylabel = "C(t)")
plot!(otoc_nW,label = label2)
plot!(OTOCapproximated, label = "exact")



#plot!([real.(Coperator(t))/D for t in 0:len])


nameplot =  @sprintf "OTOC N=%3.0i ξ = %1.2f ϵ = %1.3f length = 5e2 Coperator" N ξ ϵ

savefig(nameplot)


=#