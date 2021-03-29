using Plots
using PyCall
using Printf
using StatsBase
using LinearAlgebra


include("SpectrumFunctions.jl")
include("OTOCs.jl")


ξ = 0.001
ϵ = 0.0
N = 15

#
# pred spustenim julie 1.3 export do terminalu ~$ JULIA_NUM_THREADS=15
#

println("Num of threads: $(Threads.nthreads())") #number of threads 


len = Int64(5e4) #time


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




