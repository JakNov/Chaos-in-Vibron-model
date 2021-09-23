using Plots
using PyCall
using Printf
using LinearAlgebra
using StatsBase
using FFTW

include("SpectrumFunctions.jl")

function OTOCevaluation(W::Array{Float64,2},V::Array{Float64,2},t::Float64,spectrum::Array{Float64,1})

    Dim = length(spectrum)
    W_t = Diagonal(exp.( 1im *t*spectrum)) * W * Diagonal(exp.( -1im *t*spectrum))

    WV = W_t * V
    VW = V * W_t

    OTOC = 2/Dim * ((WV - VW)*VW - (WV - VW)*WV)
    return OTOC
end


function OTOCgeneral(W, V, Hamiltonian; len::Int64 = 5000)
    

    EigSystem = eigen(Hamiltonian)
    spectrum = EigSystem.values
    vectors = EigSystem.vectors
    
    Dim = length(spectrum)
    
    S = vectors
    W_matrix = inv(S)*W*S #n ve vlastni bazi hamiltonianu
    V_matrix = inv(S)*V*S
    

    OTOCarray = Array{Float64,2}(undef,len,Dim)
    OTOCarrayNum = Array{Float64,2}(undef,len,Dim)

    OTOCeig = Array{Float64,2}(undef,len,Dim)
    OTOCtrace = Array{Float64,1}(undef,len)


    Threads.@threads for i in 1:len
            println(i)
             
            OTOCmatrix = OTOCevaluation(W_matrix,V_matrix,(i-1)/10,spectrum)
            OTOCmatrixNum = S*OTOCmatrix*inv(S) 
            OTOCtrace[i] = real(tr(OTOCmatrix))
            #println(OTOCtrace[i])
            EigSystemOTOC = eigen(OTOCmatrix)
            spectrumOTOC = EigSystemOTOC.values
            #println(sum(spectrumOTOC))
            #println(Dim)
            for j in 1:Dim
                OTOCarray[i,j] = OTOCmatrix[j,:]'*OTOCmatrix[j,:]
                OTOCarrayNum[i,j] = OTOCmatrixNum[j,:]'*OTOCmatrixNum[j,:]

                OTOCeig[i,j] = real(spectrumOTOC[j])
            end
            
            #println(sum(OTOCeig[i,:]))
        end

    return OTOCarray, OTOCtrace, OTOCeig, OTOCarrayNum
end


ξ = 0.5
ϵ = 0.0
N = 20


H = Hamiltonian_Nnl(ξ,ϵ,N)
W2 = W2_Nnl(N)
n = n_Nnl(N)
Dx = (Dp_Nnl(N) + Dm_Nnl(N))/2   

ot1, ot2, ot3, ot4 = OTOCgeneral(n, W2, H)
sm = sum(ot1,dims = 2)
sm2 = sum(ot3,dims = 2)


pyplot()
PyPlot.pygui(true)


sg = ot1[:,30]
# Number of points 
len = length(sg) 
# Sample period
Ts = 1 / (1.0 *len) 
# Start time 
t0 = 0 
tmax = t0 + len * Ts
# time coordinate
t = t0:Ts:tmax

# signal 
signal = sg#sm#sin.(2π * 60 .* t) + 2*sin.(2π * 100 .* t)  # sin (2π f t) 

# Fourier Transform of it 
F = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

# plots 
pyplot()
PyPlot.pygui(true)
time_domain = plot(t, signal, title = "Signal")
freq_domain = plot(freqs, abs.(F))#, title = "Spectrum", xlim=(1, +500),ylim=(0,1e8))#,yaxis=:log) 
freq_domain = scatter!(freqs, abs.(F), ms=0.2 )#, title = "Spectrum", xlim=(1, +500),ylim=(0,1e8))#,yaxis=:log) 
mn = mean(abs.(F))
freq_domain = plot!(freqs,[mn for  i in 1:length(abs.(F))])

data = abs.(F)
data = data[Int64(length(data)/2) + 2:end]
data = filter( x -> x > mn, data)
f =sum(data.^4) /(data'*data)^2

plot(time_domain, freq_domain, layout = (2,1))


#plot(ot1[:,1])

