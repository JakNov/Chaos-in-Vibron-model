using Plots
using PyCall
using Statistics
using ProgressBars
using Printf

function ReadData(n::Int64,x::Float64,E0::Float64,C::Float64,Step::Float64,minn::Float64,maxx::Float64,path::String,Exeptions::Array{Float64};A::Float64 = 0.1,B::Float64 = -0.8,IsX::Bool = false,IsE::Bool = false,IsC::Bool=false)

    Surface = Float64[]
    Λ0 = Float64[]
    χS = Float64[]
    χSΛ0 = Float64[]
    Z_Surface = Float64[]
    Z_Λ0 = Float64[]
    Z_χS = Float64[]
    Z_χSΛ0 = Float64[]

    len = Int64(round((maxx-minn)/Step))

    
    for i in 0:len
        if IsX == true
            zi = x + Step*i
        elseif IsE == true
            zi = E0 + Step*i  
        elseif IsC == true
            zi = C + Step*i
        end  

        #println(Exeptions)
        #println(any(z->abs(z) <= 1e-5, Exeptions.-zi))
        if !any(z->abs(z) <= 1e-5, Exeptions.-zi) || length(Exeptions) == 0
            if IsX == true
                #zi = x + Step*i
                name = @sprintf "Results_n=%1.0i,x=%1.2f,E0=%1.2f,C=%1.2f,B=%1.1f,A=%1.1f" n zi E0 C B A
            elseif IsE == true
                #zi = E0 + Step*i
                name = @sprintf "Results_n=%1.0i,x=%1.1f,E0=%1.2f,C=%1.2f,B=%1.1f,A=%1.1f" n x zi C B A
            elseif IsC == true
                #zi = C + Step*i
                name = @sprintf "Results_n=%1.0i,x=%1.2f,E0=%1.2f,C=%1.2f,B=%1.1f,A=%1.1f" n x E0 zi B A
            else
                println("You have done it yourself!")
            end

            name = path*name*".txt" 
            for line in eachline(name)     
                line = replace(line, "[" => "")   
                line = replace(line, "]" => "")
                line = replace(line, "," => "")
                elements = split(line)
                    
                chis = (parse(Float64, elements[3]))
                s =(parse(Float64, elements[1]))
                l0 = (parse(Float64, elements[2]))
                chisl0 = (parse(Float64, elements[4]))
                
                if !isnan(s)     
                    append!(Surface, s)
                    append!(Z_Surface, zi)
                end
    
                if !isnan(l0)    
                    append!(Λ0, l0)
                    append!(Z_Λ0, zi)
                end
                    
                if !isnan(chis) && abs(chis) < 1.0    
                    append!(χS, chis)
                    append!(Z_χS,zi)
                end
        
                if !isnan(chisl0) && abs(chisl0) < 1.0  
                    append!(χSΛ0, chisl0)
                    append!(Z_χSΛ0, zi)
                elseif s>0  && abs(chis) < 1.0 && !isnan(l0) 
                        
                    append!(χSΛ0, chis/1.0)
                    append!(Z_χSΛ0, zi)
                end
            end
        end
    end
          
    return   Surface,Λ0,χS,χSΛ0,Z_Surface,Z_Λ0,Z_χS,Z_χSΛ0
         
    
end


function Read_File(n::Int64,path::String)

    Data = [[0.0] for i in 1:n, j in 1:n, k in 1:n]

    open(path) do file

    for i in 1:n #py
        for j in 1:n #px
            for k in 1:n #y

                line = readline(file)
                line = replace(line, "[" => "")   
                line = replace(line, "]" => "")
                line = replace(line, "," => "")
                elements = split(line)

                Data[k,j,i] = [parse(Float64, element) for element in elements]

 
            end
        end
    end
    end #close file

    return Data
end


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