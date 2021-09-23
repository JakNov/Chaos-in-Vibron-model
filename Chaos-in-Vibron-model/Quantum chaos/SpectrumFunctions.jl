using Plots
using PyCall
using LinearAlgebra
using SparseArrays
#using ArnoldiMethod
#using Arpack

function Commutator(W,V)
    return W*V - V*W
end

function AntiCommutator(W,V)
    return W*V + V*W
end


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
 #  Basis
 #
 function Basis_Nnl(N::Int64)
    B1 = []
    #D = Int((N + 1)*(N +2)/2) 
    for n in 0:N
        max = (n + 1)
        for i in 1:max
            append!(B1,[N,n,-n + (i-1)*2])
        end
    end

    B1 = reshape(B1,(3,:))

    return B1

end

function Basis_Nln(N::Int64)
    B2 = []
    #D = Int((N + 1)*(N +2)/2) 
    for j in 1:2N+1
        l = -N + j -1
        n = abs(l)
        while n <= N
            append!(B2,[N,l,n])
            n = n + 2
        end
    end

    
    B2 = reshape(B2,(3,:))
    
    return B2
end

function BasisList_Nnl(N::Int64)
    B1 = []
    #D = Int((N + 1)*(N +2)/2) 
    for n in 0:N
        max = (n + 1)
        for i in 1:max
            push!(B1,[N,n,-n + (i-1)*2])
        end
    end

    #B1 = reshape(B1,(3,:))

    return B1

end

function BasisList_Nln(N::Int64)
    B2 = []
    #D = Int((N + 1)*(N +2)/2) 
    for j in 1:2N+1
        l = -N + j -1
        n = abs(l)
        while n <= N
            push!(B2,[N,n,l])
            n = n + 2
        end
    end

    
    #B2 = reshape(B2,(3,:))
    
    return B2
end

function BaisiListlml(N::Int64) 

    B1 =  BasisList_Nln(N)
    Blml = []
       

    D = Int((N + 1)*(N +2)/2)

    num = 0.0
    int = 0.0

    if N%2 ==0
        global num = D - (N/2 + 1)
        global int = (N/2 + 1)
    else
        global num = D - ((N-1)/2 + 1)
        global int = ((N-1)/2 + 1)
    end


    # plus prostor

    for i in 1:Int(num/2 + int)

       push!(Blml, [N,B1[i][2],-B1[i][3]]) # N n l v bazi Nln

    end

    #minus prostor
    for i in Int(num/2 + int):D-1
        push!(Blml, [N,B1[1+i-Int(num/2 + int)][2],B1[1+i-Int(num/2 + int)][3]]) # N n l v bazi Nln
    end

    return Blml

end


function TrMalmlTolmlReordered(N::Int64)

    #B1 =  BasisList_Nln(N)
    Blml = BaisiListlml(N)

    D = Int((N + 1)*(N +2)/2)

    TransitionMatrix = zeros(Float64,D,D)

    num = 0.0 #length of first subspace
    int = 0.0

    if N%2 ==0
        global num = D - (N/2 + 1)
        global int = (N/2 + 1)
    else
        global num = D - ((N-1)/2 + 1)
        global int = ((N-1)/2 + 1)
    end



    for i in 1:Int(num/2)
        if abs(Blml[i][2])%2 == 1
            ind = findfirst(x -> x == [N,Blml[i][2],-Blml[i][3]],Blml) # N n l v bazi Nln
            TransitionMatrix[ind,i] = 1
            TransitionMatrix[i,ind] = 1

        else 
            TransitionMatrix[i,i] = 1
        end
    end

    for i in Int(num/2)+1:D
        if abs(Blml[i][2])%2 == 0
            TransitionMatrix[i,i] = 1
        end

    end


    return TransitionMatrix

end

function TrMaNnlToNln(N::Int64)

    B1 =  BasisList_Nnl(N)
    B2 =  BasisList_Nln(N)
       

    D = Int((N + 1)*(N +2)/2)

    TransitionMatrix = zeros(Float64,D,D)

    for i in 1:D
        ind = findfirst(x -> x == B1[i],B2)
        TransitionMatrix[i,ind] = 1.0
    end

    return TransitionMatrix

end


function TrMaNlnTolml(N::Int64)

    B1 =  BasisList_Nln(N)
       

    D = Int((N + 1)*(N +2)/2)

    TransitionMatrix = zeros(Float64,D,D)

    num = 0.0
    int = 0.0

    if N%2 ==0
        global num = D - (N/2 + 1)
        global int = (N/2 + 1)
    else
        global num = D - ((N-1)/2 + 1)
        global int = ((N-1)/2 + 1)
    end



    for i in 1:Int(num/2)

        ind = findfirst(x -> x == [N,B1[i][2],-B1[i][3]],B1) # N n l v bazi Nln

       # println(ind)
        TransitionMatrix[i,i] = 1.0/sqrt(2)

        TransitionMatrix[ind,i] = 1.0/sqrt(2)

        TransitionMatrix[i,D-Int(num/2) + i] = -1.0/sqrt(2)
        TransitionMatrix[ind,D-Int(num/2) + i] = 1.0/sqrt(2)
    end

    for i in Int(num/2 + 1): Int(num/2 + int)
        TransitionMatrix[i,i] = 1.0
    end


    return TransitionMatrix

end


#
#   Operators in different basis
#

function n_Nnl(N::Int64)
    Basis = Basis_Nnl(N)
    n = convert(Array{Float64,1}, Basis[2,:])

    return Diagonal(n)
end

function n_Nln(N::Int64)
    Basis = Basis_Nln(N)
    n = convert(Array{Float64,1}, Basis[3,:])

    return Diagonal(n)
end

function l_Nnl(N::Int64)
    Basis = Basis_Nnl(N)
    l = convert(Array{Float64,1}, Basis[3,:])

    return Diagonal(l)
end

function l_Nln(N::Int64)
    Basis = Basis_Nln(N)
    l = convert(Array{Float64,1}, Basis[2,:])

    return Diagonal(l)
end


function ns_Nnl(N::Int64)
    Basis = Basis_Nnl(N)
    ns = convert(Array{Float64,1}, Basis[1,:].-Basis[2,:])

    return Diagonal(ns)
end

function W2_Nnl(N::Int64)
    # |N n l> basis

    Basis = Basis_Nnl(N)

    D = Int((N + 1)*(N +2)/2) 

    W2 = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  
        W2[i,i]=   1/2 * (P2(N,n+1,l-1)*M1(N,n,l) + P1(N,n-1,l-1)*M2(N,n,l) + M2(N,n+1,l+1)*P1(N,n,l) + M1(N,n-1,l+1)P2(N,n,l) +2*l^2)
    
        if n + 2 <= N 
            W2[i,i+2*n+4]=  1/2 * (P1(N,n+1,l-1)*M1(N,n,l) + M1(N,n+1,l+1)*P1(N,n,l))
            W2[i+2*n+4,i] = W2[i,i+2*n+4]
        end
    end

    return W2
end

#
#   všechno je transponované, takže x -> y - možná už ne tohle bude chtít kontrolu
#
#=
function Dx_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Dx = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+2) <= D
            Dx[i,i+ (n+2)]= P1(N,n,l)/2
            Dx[i+ (n+2),i] = Dx[i,i+ (n+2)]
        end
    
        if i + (n+1) <= D
            Dx[i,i+ (n+1)] = M1(N,n,l)/2
            Dx[i+ (n+1),i] = Dx[i,i+ (n+1)]
        end
    
    end

    return Dx
end=#

function Dp_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Dp = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+2) <= D 
            Dp[i,i+ (n+2)]= P1(N,n,l)
        end
        
        if i - n >0
            Dp[i,i - n] = P2(N,n,l)
        end
    end

    return Dp
end

function Qp_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Qp = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+2) <= D 
            Qp[i,i+ (n+2)]= 1/2 * sqrt(n+l+2)
            #println([[n,l],[Basis[2,i+n+1],Basis[3,i+n+2]]])
        end
        
        if i - n >0
            Qp[i,i - n] = -1/2 * sqrt(n-l)
            #println([[n,l],[Basis[2,i-n],Basis[3,i-n]]])

        end
    end

    return Qp
end

function Qm_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Qm = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+1) <= D 
            Qm[i,i+ (n+1)]= 1/2 * sqrt(n-l+2)
            #println([[n,l],[Basis[2,i+n+1],Basis[3,i+n+2]]])
        end
        
        if i - n >1
            Qm[i,i - n-1] = -1/2 * sqrt(n+l)
            #println([[n,l],[Basis[2,i-n],Basis[3,i-n]]])

        end
    end

    return Qm
end

function Pp_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Pp = zeros(ComplexF64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+1) <= D 
            Pp[i,i+ (n+1)]= -im/2 * sqrt(n-l+2)
            #println([[n,l],[Basis[2,i+n+1],Basis[3,i+n+2]]])
        end
        
        if i - n >1
            Pp[i,i - n-1] = -im/2 * sqrt(n+l)
            #println([[n,l],[Basis[2,i-n],Basis[3,i-n]]])

        end
    end

    return Pp
end

function Pm_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Pm = zeros(ComplexF64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+2) <= D 
            Pm[i,i+ (n+2)]= -im/2 * sqrt(n+l+2)
            #println([[n,l],[Basis[2,i+n+1],Basis[3,i+n+2]]])
        end
        
        if i - n >0
            Pm[i,i - n] = -im/2 * sqrt(n-l)
            #println([[n,l],[Basis[2,i-n],Basis[3,i-n]]])

        end
    end

    return Pm
end

function Rp_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Rp = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+2) <= D 
            Rp[i,i+ (n+2)]= P1(N,n,l)
        end
        
        if i - n >0
            Rp[i,i - n] = -P2(N,n,l)
        end
    end

    return Rp
end

function Dm_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Dm = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+1) <= D 
            Dm[i,i+ (n+1)]= M1(N,n,l)
        end
    
        if i - n > 1
            Dm[i,i- n-1] = M2(N,n,l)
        end
    
    end

    return Dm
end

function Rm_Nnl(N::Int64)
    # |N n l> basis
    Basis = Basis_Nnl(N)
  

    D = Int((N + 1)*(N +2)/2) 

    Rm = zeros(Float64,D,D)
    for i in 1:D
        n,l = Int64(Basis[2,i]),Int64(Basis[3,i])  

        if i + (n+1) <= D 
            Rm[i,i+ (n+1)]= -M1(N,n,l)
        end
    
        if i - n > 1
            Rm[i,i- n-1] = M2(N,n,l)
        end
    
    end

    return Rm
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

function Dx_Nnl(ξ::Float64,ϵ::Float64,N::Int64)
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
    
        if i + (n+2) <= D
            H_1[i,i+ (n+2)]= -ϵ/2 *P1(N,n,l)
            H_1[i+ (n+2),i] = H_1[i,i+ (n+2)]
        end
    
        if i + (n+1) <= D
            H_1[i,i+ (n+1)] = -ϵ/2 *M1(N,n,l)
            H_1[i+ (n+1),i] = H_1[i,i+ (n+1)]
        end
    
    end

    return H_1
end

function Dx_Nln(ξ::Float64,ϵ::Float64,N::Int64)
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
                H_2[i,i+num_l] = +ϵ/2 *P2(N,n,l)
                H_2[i+num_l,i] = H_2[i,i+num_l]
            elseif  n != abs(l)
                H_2[i,i+num_l-1] = +ϵ/2 *P2(N,n,l)
                H_2[i+num_l-1,i] = H_2[i,i+num_l-1]
            end
        end
    
    end
    
    return H_2
end


function HamiltonianR_Nln(ξ::Float64,ϵ::Float64,N::Int64)
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
            H_2[i,i+1] = +1/2 * ξ/(N-1) * (P1(N,n+1,l-1)*M1(N,n,l) + M1(N,n+1,l+1)*P1(N,n,l))
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
                H_2[i,i+num_l] = +ϵ/2 *P2(N,n,l)
                H_2[i+num_l,i] = H_2[i,i+num_l]
            elseif  n != abs(l)
                H_2[i,i+num_l-1] = +ϵ/2 *P2(N,n,l)
                H_2[i+num_l-1,i] = H_2[i,i+num_l-1]
            end
        end
    
    end
    
    return H_2
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

#
#
#   Subspaces
#
#

function HdNlnSubspaces(N::Int64)

    U1 = TrMaNlnTolml(N)
    U2 = TrMalmlTolmlReordered(N)

    TrMatrix = U1*U2
   # Vectors = inv(TrMatrix)

    D = Int((N + 1)*(N +2)/2)

    num = 0.0 #length of first subspace
    int = 0.0 

    if N%2 ==0
        global num = D - (N/2 + 1)
        global int = (N/2 + 1)
    else
        global num = D - ((N-1)/2 + 1)
        global int = ((N-1)/2 + 1)
    end

    Subspace1 = TrMatrix[1:end,1:Int(num/2 + int)]   #plus prostor
    Subspace2 = TrMatrix[1:end,Int(num/2 + int)+1:end]   #minus prostor



    return Subspace1,Subspace2


end

function VectorProjectionsNln(vec,N::Int64)

    Subspace1,Subspace2 = HdNlnSubspaces(N)
    c1,c2 = vec'*Subspace1,vec'*Subspace2

    return c1*c1',c2*c2'
end
