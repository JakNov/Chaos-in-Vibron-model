using Plots
using PyCall
using LinearAlgebra
using SparseArrays
#using ArnoldiMethod
#using Arpack



#
#   Functions used to evaluate elements in Hamiltonian
#

 function P1(N::Int64,n::Int64,l::Int64)
    sqrt((n + l +2)*(N - n))
 end

 function P2(N::Int64,n::Int64,l::Int64)
    -sqrt((N - n +1)*(n - l))
 end

 function M1(N::Int64,n::Int64,l::Int64)
    -sqrt((n - l +2)*(N - n))
 end

 function M2(N::Int64,n::Int64,l::Int64)
     sqrt((N - n +1)*(n + l))
 end

#
#   Hamiltonians in different basis
#

function Hamiltonian_Nnl(ξ::Float64,ϵ::Float64,N::Int64)
    # |N n l> basis
    B1 = []
    for n in 0:N
        max = (n + 1)
        for i in 1:max
            push!(B1,[n,-n + (i-1)*2])
        end
    end

    D = Int((N + 1)*(N +2)/2) 

    H_1 = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(B1[i][1]),Int64(B1[i][2])  
        H_1[i,i]= (1 - ξ)*n - 1/2 * ξ/(N-1)*(P2(N,n+1,l-1)*M1(N,n,l) + P1(N,n-1,l-1)*M2(N,n,l) + M2(N,n+1,l+1)*P1(N,n,l) + M1(N,n-1,l+1)P2(N,n,l) +2*l^2)
    
        if n + 2 <= N 
            H_1[i,i+2*n+4]= -1/2 * ξ/(N-1) * (P1(N,n+1,l-1)*M1(N,n,l) + M1(N,n+1,l+1)*P1(N,n,l))
            H_1[i+2*n+4,i] = H_1[i,i+2*n+4]
        end
    
        if i + (n+2) <= D
            H_1[i,i+ (n+2)]= -ϵ/2 *P1(N,n,l)
            H_1[i+ (n+2),i] = H_1[i,i+ (n+2)]
        end
    
        if i + (n+1) <= D
            H_1[i,i+ (n+1)] = -ϵ/2 *M1(N,n,l)
            H_1[i+ (n+1),i] = H_1[i,i+ (n+1)]
        end
    
    end

    return Symmetric(H_1)
end

function Hamiltonian_Nln(ξ::Float64,ϵ::Float64,N::Int64)
    # |N l n> basis

    B2 = []
    for j in 1:2N+1
        l = -N + j -1
        n = abs(l)
        while n <= N
            push!(B2,[l,n])
            n = n + 2
        end
    end

    D = Int((N + 1)*(N +2)/2) 

    H_2 = zeros(Float64,D,D)

    for i in 1:D
        l,n = Int64(B2[i][1]),Int64(B2[i][2])   
        H_2[i,i] = (1 - ξ)*n - 1/2 * ξ/(N-1)*(P2(N,n+1,l-1)*M1(N,n,l) + P1(N,n-1,l-1)*M2(N,n,l) + M2(N,n+1,l+1)*P1(N,n,l) + M1(N,n-1,l+1)P2(N,n,l) +2*l^2) 
        
    
        if n+2 <=N
            H_2[i,i+1] = -1/2 * ξ/(N-1) * (P1(N,n+1,l-1)*M1(N,n,l) + M1(N,n+1,l+1)*P1(N,n,l))
            H_2[i+1,i] = H_2[i,i+1]
        end
    
        num_l = Int64(ceil((N - abs(l) +1)/2))
        num_ln = Int64(ceil((N - abs(l-1) +1)/2))
    
        
        if i+ num_l <= D && n < N
            if l <0
                H_2[i,i+num_l+1] = -ϵ/2 *P1(N,n,l)
                H_2[i+num_l+1,i] = H_2[i,i+num_l+1]
            else
                H_2[i,i+num_l] = -ϵ/2 *P1(N,n,l)
                H_2[i+num_l,i] = H_2[i,i+num_l]
            end
        end
    
        if i+ num_l-1 <= D
            if l <0
                H_2[i,i+num_l] = -ϵ/2 *P2(N,n,l)
                H_2[i+num_l,i] = H_2[i,i+num_l]
            elseif  n != abs(l)
                H_2[i,i+num_l-1] = -ϵ/2 *P2(N,n,l)
                H_2[i+num_l-1,i] = H_2[i,i+num_l-1]
            end
        end
    
    end
    
    return Symmetric(H_2)
end

function Sparse_Nnl(ξ::Float64,ϵ::Float64,N::Int64) 
    # |N n l> basis

    D = Int((N + 1)*(N +2)/2)

    I1 = Float64[]
    J1 = Float64[]
    V1 = Float64[]


    B1 = []
    for n in 0:N
        max = (n + 1)
        for i in 1:max
            push!(B1,[n,-n + (i-1)*2])
        end
    end

    for i in 1:D
        n,l = Int64(B1[i][1]),Int64(B1[i][2])  
        append!(I1,i)
        append!(J1,i)
        append!(V1,(1 - ξ)*n - 1/2 * ξ/(N-1)*(P2(N,n+1,l-1)*M1(N,n,l) + P1(N,n-1,l-1)*M2(N,n,l) + M2(N,n+1,l+1)*P1(N,n,l) + M1(N,n-1,l+1)P2(N,n,l) +2*l^2))

        if n + 2 <= N 
            append!(I1,i)
            append!(J1,i+2*n+4)
            val = -1/2 * ξ/(N-1) * (P1(N,n+1,l-1)*M1(N,n,l) + M1(N,n+1,l+1)*P1(N,n,l))
            append!(V1,val)

            append!(I1,i+2*n+4)
            append!(J1,i)
            append!(V1,val)
        end

        if i + (n+2) <= D
            append!(I1,i)
            append!(J1,i+ (n+2))
            val = -ϵ/2 *P1(N,n,l)
            append!(V1,val)

            append!(I1,i+ (n+2))
            append!(J1,i)
            append!(V1,val)
        end
    
        if i + (n+1) <= D
            append!(I1,i)
            append!(J1,i+ n+1)
            val = -ϵ/2 *M1(N,n,l)
            append!(V1,val)

            append!(I1,i+ n+1)
            append!(J1,i)
            append!(V1,val)
        end

    end

    return sparse(I1,J1,V1)
end

function Sparse_Nln(ξ::Float64,ϵ::Float64,N::Int64) 
    # |N l n> basis

    D = Int((N + 1)*(N +2)/2)

    I2 = Float64[]
    J2 = Float64[]
    V2 = Float64[]

    B2 = []
    for j in 1:2N+1
        l = -N + j -1
        n = abs(l)
        while n <= N
            push!(B2,[l,n])
            n = n + 2
        end
    end


    for i in 1:D
        l,n = Int64(B2[i][1]),Int64(B2[i][2])   
        append!(I2,i)
        append!(J2,i)
        append!(V2,(1 - ξ)*n - 1/2 * ξ/(N-1)*(P2(N,n+1,l-1)*M1(N,n,l) + P1(N,n-1,l-1)*M2(N,n,l) + M2(N,n+1,l+1)*P1(N,n,l) + M1(N,n-1,l+1)P2(N,n,l) +2*l^2))
        
    
        if n+2 <=N
            append!(I2,i)
            append!(J2,i+1)
            val = -1/2 * ξ/(N-1) * (P1(N,n+1,l-1)*M1(N,n,l) + M1(N,n+1,l+1)*P1(N,n,l))
            append!(V2,val)
     
            append!(I2,i+1)
            append!(J2,i)
            append!(V2,val)
        end

        num_l = Int64(ceil((N - abs(l) +1)/2))
        num_ln = Int64(ceil((N - abs(l-1) +1)/2))
    
        
        if i+ num_l <= D && n < N
            if l <0
    
                append!(I2,i)
                append!(J2,i+num_l+1)
                val = -ϵ/2 *P1(N,n,l)
                append!(V2,val)
     
                append!(I2,i+num_l+1)
                append!(J2,i)
                append!(V2,val)
            else
                append!(I2,i)
                append!(J2,i+num_l)
                val = -ϵ/2 *P1(N,n,l)
                append!(V2,val)
     
                append!(I2,i+num_l)
                append!(J2,i)
                append!(V2,val)
            end
        end
    
        if i+ num_l-1 <= D
            if l <0
                append!(I2,i)
                append!(J2,i+num_l)
                val = -ϵ/2 *P2(N,n,l)
                append!(V2,val)
     
                append!(I2,i+num_l)
                append!(J2,i)
                append!(V2,val)
            elseif  n != abs(l)
                append!(I2,i)
                append!(J2,i+num_l-1)
                val = -ϵ/2 *P2(N,n,l)
                append!(V2,val)
     
                append!(I2,i+num_l-1)
                append!(J2,i)
                append!(V2,val)
            end
        end

    end

    return sparse(I2,J2,V2)
end