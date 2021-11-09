using Printf
using Plots
using PyCall

function Read_Results(path::String)

    Data = []
    open(path) do file 
        line = readline(file)
        line = replace(line, "[" => "")   
        line = replace(line, "]" => "")
        line = replace(line, "," => "")
        elements = split(line)
        Data = [parse(Float64, element) for element in elements]
    end


    return Data
end

function PlotData!(N::Int64,ξ::Float64,ϵ::Float64;OTOCname::String = "[qx(t),px]",path::String = "")

    pars = @sprintf "_N=%1.0i,chi=%1.2f,eps=%1.2f" N ξ ϵ
    Spectrum = Read_Results(path*"Spectrum"*OTOCname*pars*".txt")
    Mean = Read_Results(path*"Mean"*OTOCname*pars*".txt")
    Var = Read_Results(path*"Var"*OTOCname*pars*".txt")

    #pyplot()
    plotly()
    scatter(Spectrum,Mean,label = "Mean", title = "Mean"*pars)
    scatter!(Spectrum,Var,label = "Var", title = "Mean"*pars)

    namefig = @sprintf("Mean %s chi=%1.2f eps=%1.2f N=%d .html", OTOCname, ξ, ϵ, N)
    #namefig = @sprintf("Mean %s chi=%1.2f eps=%1.2f N=%d.png", OTOCname, ξ, ϵ, N)
    
    #savefig("VarMean $name chi=$ξ eps=$ϵ N=$N")
    savefig(namefig)

    VarMean = Var./Mean
    scatter(Spectrum,VarMean, title = "Var/Mean"*pars)
    #namefig = @sprintf("VarMean %s chi=%1.2f eps=%1.2f N=%d.png", OTOCname, ξ, ϵ, N)
    
    namefig = @sprintf("VarMean %s chi=%1.2f eps=%1.2f N=%d .html", OTOCname, ξ, ϵ, N)
    #savefig("VarMean $name chi=$ξ eps=$ϵ N=$N")
    savefig(namefig)

    pyplot()
    scatter(Spectrum,Mean,label = "Mean", title = "Mean"*pars)
    scatter!(Spectrum,Var,label = "Var", title = "Mean"*pars)

    namefig = @sprintf("Mean %s chi=%1.2f eps=%1.2f N=%d.png", OTOCname, ξ, ϵ, N)
    
    #savefig("VarMean $name chi=$ξ eps=$ϵ N=$N")
    savefig(namefig)

    VarMean = Var./Mean
    scatter(Spectrum,VarMean, title = "Var/Mean"*pars)
    namefig = @sprintf("VarMean %s chi=%1.2f eps=%1.2f N=%d.png", OTOCname, ξ, ϵ, N)

    savefig(namefig)
end