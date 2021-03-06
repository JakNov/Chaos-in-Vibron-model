using Distributed
workers = 2 #cores of PC

#nprocs() pocet procesu - > chceme pracovat na 4 jadrech - > chceme 4+1 procesu
if nprocs() <= workers
    addprocs(workers + 1 - nprocs())
end

@everywhere include("/home/jakub/Julia/Bc/remote/PS2.jl")


   
@everywhere function MakePoincare(n::Int64,x::Float64,E0::Float64,C::Float64;path::String = "/home/jakub/Julia/Bc/PoincareVibron_",A::Float64 = 0.1,B::Float64 = -0.8 )
    
    results,Points,Lyapunovs = PoincareSection(n,x,E0,C,A,B)
    open(path*"Lya_n=$n,x=$x,E0=$E0,C=$C,B=$B,A=$A"*".txt", "a") do lya   
        for i in 1:n #py
            for j in 1:n #px
                for k in 1:n #y
                    println(lya,Lyapunovs[k,j,i])
                end
            end
        end
    end

    open(path*"Points_n=$n,x=$x,E0=$E0,C=$C,B=$B,A=$A"*".txt", "a") do poi        
        for i in 1:n #py
            for j in 1:n #px
                for k in 1:n #y
                    println(poi,Points[k,j,i])
                end
            end
        end
    end

    open(path*"Results_n=$n,x=$x,E0=$E0,C=$C,B=$B,A=$A"*".txt", "a") do res        
        
        println(res,results)
      
    end


    return
end


function RunMapping(n::Int64,x::Float64,E0::Float64,C_min::Float64)
   
    path= "/home/jakub/Julia/Bc/PoincareVibron_"
  

    input = [(n,x,E0,C_min + i*0.01) for i in 0:10]
    pmap((args)->MakePoincare(args...;path=path), input)    
    return
end



function ReadFiles(a::Float64;path::String = "/home/jakub/Julia/GayMum_")
    file = path*"y$a"*".txt"

    list = Array{Float64}[]
    for line in eachline(file)
        line = replace(line, "[" => "")
        line = replace(line, "]" => "")
        line = replace(line, "," => "")
        elements = split(line)
        append!(list,[[parse(Float64, elements[1]), parse(Float64, elements[2]), parse(Float64, elements[3]), parse(Float64, elements[4])]])
    end

    return list
 
end