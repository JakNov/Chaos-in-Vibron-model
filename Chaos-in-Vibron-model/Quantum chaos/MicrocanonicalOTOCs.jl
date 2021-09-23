using Plots
using PyCall
using Printf
using LinearAlgebra
using StatsBase
using FFTW
using Statistics

include("SpectrumFunctions.jl")

function OTOCevaluation(W,V,t::Float64,spectrum::Array{Float64,1})

    Dim = length(spectrum)
    W_t = Diagonal(exp.( 1im *t*spectrum)) * W * Diagonal(exp.( -1im *t*spectrum))

    WV = W_t * V
    VW = V * W_t

    OTOC = 2/Dim * ((WV - VW)*VW - (WV - VW)*WV)
    return OTOC
end


function OTOCgeneral(W, V, Hamiltonian; len::Int64 = 10000)
    

    EigSystem = eigen(Hamiltonian)
    spectrum = EigSystem.values
    vectors = EigSystem.vectors
    
    Dim = length(spectrum)
    
    S = vectors
    W_matrix = inv(S)*W*S #n ve vlastni bazi hamiltonianu
    V_matrix = inv(S)*V*S
    

    OTOCarrayEnergyScale = Array{Float64,2}(undef,len,Dim)
    OTOCarrayParticleScale = Array{Float64,2}(undef,len,Dim)


    Threads.@threads for i in 1:len
            println("$i" )
             
            OTOCmatrix = OTOCevaluation(W_matrix,V_matrix,(i-1)/10,spectrum)
            OTOCmatrixP = S*OTOCmatrix*inv(S)


            for j in 1:Dim
                OTOCarrayEnergyScale[i,j] = OTOCmatrix[j,:]'*OTOCmatrix[j,:]
                OTOCarrayParticleScale[i,j] = OTOCmatrixP[j,:]'*OTOCmatrixP[j,:]

            end
            
            #println(sum(OTOCeig[i,:]))
        end

    return OTOCarrayEnergyScale, OTOCarrayParticleScale
end

function OscilationValue(signal::Array{Float64,1},name::String)
    len = length(signal) 
    # Sample period
    Ts = 1 / (1.0 *len) 
    # Start time 
    t0 = 0 
    tmax = t0 + len * Ts
    # time coordinate
    t = t0:Ts:tmax
    
    # signal 
    #signal = signal#sm#sin.(2π * 60 .* t) + 2*sin.(2π * 100 .* t)  # sin (2π f t) 
    
    # Fourier Transform of it 
    F = fft(signal) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    
    mn = mean(abs.(F))
    data = abs.(F)
    data = data[Int64(length(data)/2) + 10:end]
    pyplot()
    #PyPlot.pygui(false)
    
    #println(data[1:10])
    #data = filter( x -> x > mn, data)
    #println(data[1:10])
    
    
    dat = copy(data)

    for i in 2:length(data)-1
        if data[i-1] >= data[i] || data[i] <= data[i+1]
            dat[i]=0
        end
              
    end

    data = data/sqrt((data'*data))
    dat = dat/sqrt((dat'*dat))
   
    f =sum(dat.^4)


    #
    #   plot jednotlivých hladin
    #

    #=p1 = plot(signal, title = name)
    p2 = plot(dat, title = "$f")
    plot(p1,p2,layout = (2,1))

    savefig(name)=#
    return log(1/f)
end


function SetOTOC(N::Int64,ξ::Float64,ϵ::Float64)

W2 = W2_Nnl(N)
n = n_Nnl(N)
l = l_Nnl(N)
Dx = (Dp_Nnl(N) + Dm_Nnl(N))/2   
np = (n + l)/2
nm = (n - l)/2

Qp = Qp_Nnl(N)
Qm = Qp_Nnl(N)

Pp = Pp_Nnl(N)
Pm = Pp_Nnl(N)



pyplot()

H = Hamiltonian_Nnl(ξ,ϵ,N)
#otE, otP = OTOCgeneral(n, W2, H)
#otE, otP = OTOCgeneral(np, nm, H)
otE, otP = OTOCgeneral(Qp, Pp, H)

#otE, otP = OTOCgeneral(np, np, H)


dim = Int64((N+1)*(N+2)/2)
makingpictures = [OscilationValue(otE[:,i],"EnergyScale chi=$ξ eps=$ϵ N=$N $i") for i in 1:dim]
plot(makingpictures)
savefig("EnergyScale Qp(t)Pp chi=$ξ eps=$ϵ N=$N")

MeanScale = [mean(otE[end - Int64(1e3) : end,i]) for i in 1:dim]
VarScale = [sqrt(var(otE[end - Int64(1e3) : end,i])) for i in 1:dim]

plotly()
scatter(MeanScale,yerror=VarScale)
savefig("MeanScaleScatter Qp(t)Pp chi=$ξ eps=$ϵ N=$N")

pyplot()
smooth = [sum(makingpictures[i:i+10])/10 for i in 1:dim-10]
plot(smooth)
savefig("EnergyScale Qp(t)Pp smooth10 chi=$ξ eps=$ϵ N=$N")

makingpictures = [OscilationValue(otP[:,i],"ParticleScale chi=$ξ eps=$ϵ N=$N $i") for i in 1:dim]
plot(makingpictures)
savefig("ParticleScale Qp(t)Pp chi=$ξ eps=$ϵ N=$N")

BasisNnl = BasisList_Nnl(N)
Particles = zeros(Float64, N+1, N+1)
for i in 1:dim
    vec = BasisNnl[i]
    println(vec)
    np = (vec[2]+vec[3])/2 +1
    nm = (vec[2]-vec[3])/2 +1
    Particles[Int64(np),Int64(nm)] = makingpictures[i]
    println([np,nm])
end

heatmap(Particles)
savefig("ParticleScale Qp(t)Pp heatmap chi=$ξ eps=$ϵ N=$N")

MeanScale = [mean(otP[end - Int64(1e3) : end,i]) for i in 1:dim]
VarScale = [sqrt(var(otP[end - Int64(1e3) : end,i])) for i in 1:dim]

ParticlesMean = zeros(Float64, N+1, N+1)
ParticlesVar = zeros(Float64, N+1, N+1)

for i in 1:dim
    vec = BasisNnl[i]
    println(vec)
    np = (vec[2]+vec[3])/2 +1
    nm = (vec[2]-vec[3])/2 +1
    ParticlesMean[Int64(np),Int64(nm)] = MeanScale[i]
    ParticlesVar[Int64(np),Int64(nm)] = VarScale[i]
    println([np,nm])
end

heatmap(ParticlesMean)
savefig("ParticleMean Qp(t)Pp heatmap chi=$ξ eps=$ϵ N=$N")

heatmap(ParticlesVar)
savefig("ParticleVar Qp(t)Pp heatmap chi=$ξ eps=$ϵ N=$N")
end

#SetOTOC(5,0.75,0.06)

SetOTOC(30,0.75,0.0)
#SetOTOC(15,0.75,0.005)
#SetOTOC(30,0.75,0.1)
#SetOTOC(30,0.75,0.2)
#SetOTOC(30,0.75,0.3)
#SetOTOC(30,0.75,0.4)
#SetOTOC(30,0.75,0.5)
#SetOTOC(30,0.75,0.6)
#SetOTOC(30,0.75,0.7)
#SetOTOC(30,0.75,0.8)






