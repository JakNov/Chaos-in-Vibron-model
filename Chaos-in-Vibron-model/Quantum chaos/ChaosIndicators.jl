using Plots
using LinearAlgebra
using LsqFit
using SpecialFunctions
using StatsBase
using Polynomials

Brody(s,β) = (β[1]+1)*gamma((β[1]+2)/(β[1]+1))^(β[1]+1)*s.^(β[1]) .* exp.(-gamma((β[1]+2)/(β[1]+1))^(β[1]+1)*s.^(β[1]+1))
Poisson(s) = exp(-s) 
Wigner(s) = pi*0.5*s*exp(- pi *s^2/4)

function CutDegeneracy(spectrum::Array{Float64}) #
    unique!(spectrum)
    spc = Float64[]
    append!(spc,spectrum[1])

    for i in 1:length(spectrum)-1
        if abs(spectrum[i+1] - spc[end]) > 1e-13 #1e-14 errors in spectrum
            append!(spc,spectrum[i+1])
        end  
    end

    return spc
end

function NND(spectrum::Array{Float64}) # returns histogram of NND

    #spectrum = CutDegeneracy(spectrum)
    
    len = length(spectrum)
    xs = LinRange(1.0,len,len)
    f = polyfit(spectrum,xs,7)
    spc_cut = spectrum[50:end-50]
    
    unfolded = f.(spc_cut)
    unfolded_stretched = (unfolded .- unfolded[1]) / (unfolded[end] - unfolded[1]) * (length(unfolded) - 1)
    spacing = unfolded_stretched[2:end] .- unfolded_stretched[1:end-1]
    len = length(spacing)
    filter!(e ->e >= 0,spacing)
  
    #println(spacing)
    dif = len - length(spacing)
    if dif != 0
        println("Num of negative")
        println(dif)
    end


    nbins = 75
    h_spacing = fit(Histogram,spacing,nbins = nbins)
    h_spacing = normalize(h_spacing, mode=:pdf)

    return h_spacing

end

function BrodyParameter(HistogramSpacing)

    tdata = collect(HistogramSpacing.edges[1])[1:end-1]
    ydata = HistogramSpacing.weights[1:end]
    p0 = [0.5]

    fitt = curve_fit(Brody, tdata, ydata, p0)
    param = fitt.param
 
    return param[1],stderror(fitt)[1]
end


function IPR(vectors::Array{Float64,2})

    ξ1_E = []
    ξ2_E = []
    D = length(vectors[:,1])
    ξ_E_deloc = D/3
    
    for i in 1:D
        append!(ξ1_E,sum(vectors[:,i].^4)^(-1))
    end
    
    IPR_ξ = sum(ξ1_E)/(D * ξ_E_deloc)

    return IPR_ξ
end



function η(spectrum::Array{Float64})
    d = spectrum[2:end] - spectrum[1:end-1]
    min_distribution = Float64[]
    for i in 1:length(d)-1
        r = d[i]/d[i+1]
        append!(min_distribution, min(1/r,r))
    end

    I_p = 0.386
    I_w = 0.586

    η = (mean(min_distribution) - I_p) / (I_w - I_p)

    return η
end