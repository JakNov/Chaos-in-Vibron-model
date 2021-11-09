using Plots
using PyCall
using Printf
using LinearAlgebra
using StatsBase
using FFTW
using Statistics

#
#
#  
#
#

include("SpectrumFunctions.jl")

function OTOCevaluation(W,V,t::Float64,spectrum::Array{Float64,1})

    Dim = length(spectrum)
    W_t = Diagonal(exp.( 1im *t*spectrum)) * W * Diagonal(exp.( -1im *t*spectrum))

    WV = W_t * V
    VW = V * W_t

    OTOC = 2/Dim * ((WV - VW)*VW - (WV - VW)*WV)
    return OTOC
end


function OTOCVarMean(W, V, Hamiltonian; len::Int64 = 1500, Min = 1e7, Max = 1e8)
    

    EigSystem = eigen(Hamiltonian)
    spectrum = EigSystem.values
    vectors = EigSystem.vectors
    
    Dim = length(spectrum)
    
    S = vectors
    W_matrix = inv(S)*W*S #n ve vlastni bazi hamiltonianu
    V_matrix = inv(S)*V*S
    

    OTOCarrayEnergyScale = Array{Float64,2}(undef,len,Dim)
    VarMeanEnergyScale = Array{Float64,1}(undef,Dim)


    Points = sort(rand(len))
    Interval = (Max - Min)

    Threads.@threads for i in 1:len
            println("$i" )
             
            OTOCmatrix = OTOCevaluation(W_matrix,V_matrix,Min + Points[i]*Interval,spectrum)

            for j in 1:Dim
                OTOCarrayEnergyScale[i,j] = OTOCmatrix[j,:]'*OTOCmatrix[j,:]

            end
            
        end

        for j in 1:Dim
            VarMeanEnergyScale[j] = sqrt(var(OTOCarrayEnergyScale[:,j]))/mean(OTOCarrayEnergyScale[:,j])
        end

    return VarMeanEnergyScale, spectrum
end


function SetOTOCFast(N::Int64,ξ::Float64,ϵ::Float64)

    W2 = W2_Nnl(N)
    n = n_Nnl(N)
    
    l = l_Nnl(N)
    Dx = (Dp_Nnl(N) + Dm_Nnl(N))/2   
    np = (n + l)/2
    nm = (n - l)/2
    
    Qp = Qp_Nnl(N)
    Qm = Qm_Nnl(N)
    
    Pp = Pp_Nnl(N)
    Pm = Pm_Nnl(N)
    
    qa = -sqrt(2)*(Qp - Qm)
    pa = -sqrt(2)*(Pp - Pm)
    
    H = Hamiltonian_Nnl(ξ,ϵ,N)
    #VarMeanEnergy, spectrum = OTOCVarMean(n, W2, H)
    #name = "[n(t),W2]"

    VarMeanEnergy, spectrum = OTOCVarMean(np,nm,H)
    name = "[np(t),np]"

    #println(typeof(name))
    #println(name)
    
   # plotly()
    pyplot()
    scatter(spectrum,VarMeanEnergy, title = "Var / Mean $name chi=$ξ eps=$ϵ N=$N")
    namefig = @sprintf("VarMean %s chi=%1.2f eps=%1.2f N=%d", name, ξ, ϵ, N)
    #savefig("VarMean $name chi=$ξ eps=$ϵ N=$N")
    savefig(namefig)

    return spectrum, VarMeanEnergy


end
    

function ResizeSpectrum(min::Float64, max::Float64,N::Int64, Spectrum::Array{Float64,1},Indicator)
    IndicatorEquidistant = zeros(N)
    step = (max-min)/N

    len =length(Spectrum) 
    last = 1 
    for i in 1:N
        now = last
        while now < len && Spectrum[now] < min + i*step
            #println(now)
            IndicatorEquidistant[i] += Indicator[now]

            now +=1
        end

        if now != last 
            IndicatorEquidistant[i] /= (now - last)
        end

        last = now

        #println("last $last")
    end

    return IndicatorEquidistant
end

function PlotHeatmaps(ϵ::Float64)

    otoc = "[n(t),W2]"
    otoc = "[np(t),np]"


    N = 30
    #ϵ = 0.2
    len = Int64((N+2)*(N+1)/2)
    Num = 20
    Spectra = Float64[]
    Indicators = Float64[]

    for i in 0:Num

        spectrum, IndicatorE = SetOTOCFast(N,0.0 + i/(Num),ϵ)

        append!(Spectra,spectrum)
        append!(Indicators, IndicatorE)

    end

    Max = max(Spectra...)
    Min = min(Spectra...)

    IndicatorEquid = Float64[]

    SpcLen = (len - len%5)/5
    x = 0.0:1.0/Num:1.0
    y = Min:(Max-Min)/SpcLen:Max

    for i in 0:Num

        IndEquid = ResizeSpectrum(Min,Max,Int64(SpcLen),Spectra[i*len + 1 : (i+1)*len],Indicators[i*len + 1: (i+1)*len])
        append!(IndicatorEquid,IndEquid)



    end

    pyplot()
    heatmap(reshape(IndicatorEquid,(:,Num+1)))
    savefig("HeatmapSpectrum VarMean $otoc N = $N eps = $ϵ")
end

#PlotHeatmaps(0.0)
#PlotHeatmaps(0.1)
#PlotHeatmaps(0.2)
#PlotHeatmaps(0.3)
#PlotHeatmaps(0.4)
#PlotHeatmaps(0.5)
#PlotHeatmaps(0.6)




