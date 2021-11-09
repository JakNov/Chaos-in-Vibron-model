using DifferentialEquations
using LinearAlgebra
using Plots
using PyCall
using Random 
using Statistics
using Polynomials
using Roots

#   SIGMA
function Σ(x::Float64,y::Float64,px::Float64,py::Float64) #returns this useful sum
    return (x^2 + y^2 + px^2 + py^2)
end


#   ENERGY VIBRON
function E_vibron(u::Array{Float64},p::Array{Float64}) #returns energy of point in vibron phase space
    x,y,px,py = u[1:4]
    A,B,C = p[1:3]
    Σ = x^2 + y^2 + px^2 + py^2
 
    return A*Σ + B*((px^2 + py^2)*(2-Σ) + (x*py - y*px)^2) + C*py*sqrt(2-Σ)
end


#   EVALUATION OF SURFACE
function Surface_linear(P::Array{Float64},par::Array{Float64}) #returns the surface of tangent plane to surface at point P in the box given by i,j,k coordinates 
    x,A,B,C,E0,i,j,k,n = par[1:9]

    # coefficients in Taylor series for x = const. section
    h1 = 2*A*P[1] - C*P[3]*P[1]/sqrt(2 - Σ(x,P[1],P[2],P[3])) + B*(-2*(P[2]^2 + P[3]^2)*P[1] - 2*P[2]*(P[3]*x - P[2]*P[1]))
    h2 = 2*A*P[2] - C*P[2]*P[3]/sqrt(2 - Σ(x,P[1],P[2],P[3])) + B*(-2*(P[2]^2 + P[3]^2)*P[2] - 2*P[1]*(P[3]*x - P[2]*P[1]) + 2*P[2]*(2 - Σ(x,P[1],P[2],P[3])))
    h3 = 2*A*P[3] - C*P[3]*P[3]/sqrt(2 - Σ(x,P[1],P[2],P[3])) + C*sqrt(2 - Σ(x,P[1],P[2],P[3])) + B*(-2*(P[2]^2 + P[3]^2)*P[3] + 2*x*(P[3]*x - P[2]*P[1]) + 2*P[3]*(2 - Σ(x,P[1],P[2],P[3])))
    K = -h1*P[1]-h2*P[2]-h3*P[3] 
       
    bound_x = [sqrt(2)*(-1 + 2/n *(k-1)),sqrt(2)*(-1 + 2/n *(k))]
    bound_y = [sqrt(2)*(-1 + 2/n *(j-1)),sqrt(2)*(-1 + 2/n *(j))]
    bound_z = [sqrt(2)*(-1 + 2/n *(i-1)),sqrt(2)*(-1 + 2/n *(i))]

    Points = Array{Float64}[]

    Plane(x::Float64,y::Float64,z::Float64) = h1*x+h2*y+h3*z + K
    Root_x(y::Float64,z::Float64) = (-K - y*h2 - z*h3)/h1
    Root_y(x::Float64,z::Float64) = (-K - x*h1 - z*h3)/h2
    Root_z(x::Float64,y::Float64) = (-K - x*h1 - y*h2)/h3

    
    S = 0.0
    for kk in 1:2
        for jj in 1:2
            for ii in 1:2
                if Plane(bound_x[kk],bound_y[jj],bound_z[ii]) == 0
                    push!(Points,[bound_x[kk],bound_y[jj],bound_z[ii]])
                end

            end
        end
    end

    for ii in 1:2
        for jj in 1:2
            x_root = Root_x(bound_y[ii],bound_z[jj])
            if bound_x[1] < x_root < bound_x[2]
                push!(Points,[x_root,bound_y[ii],bound_z[jj]])
            end
        end
    end

    for ii in 1:2
        for jj in 1:2
            y_root = Root_y(bound_x[ii],bound_z[jj])
            if bound_y[1] < y_root < bound_y[2]
                push!(Points,[bound_x[ii],y_root,bound_z[jj]])
            end
        end
    end

    for ii in 1:2
        for jj in 1:2
            z_root = Root_z(bound_x[ii],bound_y[jj])
            if bound_z[1] < z_root < bound_z[2]
                push!(Points,[bound_x[ii],bound_y[jj],z_root])
            end
        end
    end

    Points_length = length(Points)


    if Points_length >= 4
        Sort_Points = Array{Float64}[]
        push!(Sort_Points,Points[1])
        delky =  Float64[]
        splice!(Points, 1)
        for ii in 1:Points_length-1
            distances = [sqrt(dot(Points[jj] - Sort_Points[ii],Points[jj] - Sort_Points[ii])) for jj in 1:Points_length-ii]
            push!(Sort_Points,Points[argmin(distances)])
            append!(delky,min(distances...))
            splice!(Points, argmin(distances))
        end

       
        a = delky[1]
        b = delky[2]
        c = sqrt(dot(Sort_Points[1] - Sort_Points[3],Sort_Points[1] - Sort_Points[3]))
        s = (a+b+c)/2
        S += sqrt(s*(s - a)*(s - b)*(s - c))
        S1 = sqrt(s*(s - a)*(s - b)*(s - c))

        
        for ii in 3:Points_length-1
            a = delky[ii]
            b = sqrt(dot(Sort_Points[1] - Sort_Points[ii],Sort_Points[1] - Sort_Points[ii]))
            c = sqrt(dot(Sort_Points[1] - Sort_Points[ii+1],Sort_Points[1] - Sort_Points[ii+1]))
            s = (a+b+c)/2
            S += sqrt(s*(s - a)*(s - b)*(s - c)) #Heron formula
        end

    elseif Points_length == 3
        l2 = sqrt(dot(Points[1] - Points[2],Points[1] - Points[2])) 
        l3 = sqrt(dot(Points[1] - Points[3],Points[1] - Points[3]))
        l4 = sqrt(dot(Points[2] - Points[3],Points[2] - Points[3]))

        s1 = (l2 + l3 + l4)/2
       
        S = sqrt(s1*(s1 - l2)*(s1 - l3)*(s1 - l4)) 
    
    else 
        #println("uff $Points_length")
        S = 0.0
    end
  
    #control
    #if S ==0 
        #println(S)
    #end
    return S
end

#
#   NON POLYNOMIAL DERIVATION OF INITIAL CONDITIONS
#
#   returns all roots of equation given by Hamiltonian
#

function InitialCondS_Px(x::Float64,y::Float64,py::Float64,E0::Float64,A::Float64,B::Float64,C::Float64)
    Roots = Float64[]
    Ham(px) = E_vibron([x,y,px,py],[A,B,C]) - E0
    Ham2(px) = (A*Σ(x,y,px,py) + B*((px^2 + py^2)*(2 -  Σ(x,y,px,py)) + (x*py - y*px)^2)-E0)^2 - C^2*py^2*(2 - Σ(x,y,px,py))
   
    point = find_zeros(Ham2,-sqrt(2),sqrt(2))
    
    if Array{Float64,1}[] != point
        
        for i in 1:length(point)
            if abs(imag(point[i])) <1e-5 && Σ(x,y,real(point[i]),py) <= 2 && abs(Ham(real(point[i]))) <= 1e-9
                if i == 1
                    append!(Roots, Float64(real(point[i])))
                elseif abs(real(point[i]) - real(point[i-1])) > 1e-5
                    append!(Roots, Float64(real(point[i])))
                end
                
            end
            
        end
    end

   if length(Roots) > 0
        return Roots
    else 
        return missing
    end
end


function InitialCondS_Py(x::Float64,y::Float64,px::Float64,E0::Float64,A::Float64,B::Float64,C::Float64)
    Roots = Float64[]
    Ham(py) = E_vibron([x,y,px,py],[A,B,C]) - E0
    Ham2(py) = (A*Σ(x,y,px,py) + B*((px^2 + py^2)*(2 -  Σ(x,y,px,py)) + (x*py - y*px)^2)-E0)^2 - C^2*py^2*(2 - Σ(x,y,px,py))
   
    point = find_zeros(Ham2,-sqrt(2),sqrt(2))
    
    if Array{Float64,1}[] != point
        
        for i in 1:length(point)
            if abs(imag(point[i])) <1e-5 && Σ(x,y,px,real(point[i])) <= 2 && abs(Ham(real(point[i]))) <= 1e-9
                if i == 1
                    append!(Roots, Float64(real(point[i])))
                elseif abs(real(point[i]) - real(point[i-1])) > 1e-5
                    append!(Roots, Float64(real(point[i])))
                end
                
            end
            
        end
    end

   if length(Roots) > 0
        return Roots
    else 
        return missing
    end
end

function InitialCondS_X(y::Float64,px::Float64,py::Float64,E0::Float64,A::Float64,B::Float64,C::Float64)
    Roots = Float64[]
    Ham(x) = E_vibron([x,y,px,py],[A,B,C]) - E0
    Ham2(x) = (A*Σ(x,y,px,py) + B*((px^2 + py^2)*(2 -  Σ(x,y,px,py)) + (x*py - y*px)^2)-E0)^2 - C^2*py^2*(2 - Σ(x,y,px,py))
   
    point = find_zeros(Ham2,-sqrt(2),sqrt(2))
    
    if Array{Float64,1}[] != point
        
        for i in 1:length(point)
            if abs(imag(point[i])) <1e-5 && Σ(real(point[i]),y,px,py) <= 2 && abs(Ham(real(point[i]))) <= 1e-9
                if i == 1
                    append!(Roots, Float64(real(point[i])))
                elseif abs(real(point[i]) - real(point[i-1])) > 1e-5
                    append!(Roots, Float64(real(point[i])))
                end
                
            end
            
        end
    end

   if length(Roots) > 0
        return Roots
    else 
        return missing
    end
end

function InitialCondS_Y(x::Float64,px::Float64,py::Float64,E0::Float64,A::Float64,B::Float64,C::Float64)
    Roots = Float64[]
    Ham(y) = E_vibron([x,y,px,py],[A,B,C]) - E0
    Ham2(y) = (A*Σ(x,y,px,py) + B*((px^2 + py^2)*(2 -  Σ(x,y,px,py)) + (x*py - y*px)^2)-E0)^2 - C^2*py^2*(2 - Σ(x,y,px,py))
   
    point = find_zeros(Ham2,-sqrt(2),sqrt(2))
    
    if Array{Float64,1}[] != point
        
        for i in 1:length(point)
            if abs(imag(point[i])) <1e-5 && Σ(x,real(point[i]),px,py) <= 2 && abs(Ham(real(point[i]))) <= 1e-9
                if i == 1
                    append!(Roots, Float64(real(point[i])))
                elseif abs(real(point[i]) - real(point[i-1])) > 1e-5
                    append!(Roots, Float64(real(point[i])))
                end
                
            end
            
        end
    end

   if length(Roots) > 0
        return Roots
    else 
        return missing
    end
end


#   FIND 
function find(points::Array{Float64},a::Float64,b::Float64) #looks for point between a and b
    point =missing
    i = 1
    len = length(points)
    while point === missing && i <= len

        if points[i] >= a && points[i] <= b
            point = points[i]
        end
        i+=1
    end

    return point
end

# Projection of one point from the box given by i,j,k coordinates to the furface given by the vibron Hamiltonian
function Point_of_Intersection(Points::Array{Array{Float64}},par::Array{Float64})
    x,A,B,C,E0,i,j,k,n = par[1:9]

    X = Points[1]
    Y = Points[2]
    Z = Points[3]

    Intersection = missing

    if X != Y && X != Z && Y !=X 
    
        T = (X + Y + Z)/3
    
        Nvec=[(X[2]) - (Y[2])*(X[3] - Z[3]) - (X[3]) - (Y[3])*(X[2] - Z[2]); #normal vector to the plane given by X,Y,Z
        (X[3]) - (Y[3])*(X[1] - Z[1]) - (X[1]) - (Y[1])*(X[3] - Z[3]);
        (X[1]) - (Y[1])*(X[2] - Z[2]) - (X[2]) - (Y[2])*(X[1] - Z[1])]

        vec(t) = T + t*Nvec

        #urceni hranic t vazane na X,Y,Z
        minima = Float64[]
        maxima = Float64[]

        if Nvec[1] != 0
            Bounds_y = ([sqrt(2)*(-1 + 2/n *(k-1)),sqrt(2)*(-1 + 2/n *k)] - T[1]*[1.0,1.0])/Nvec[1]

            if Bounds_y[1] >= Bounds_y[2]
                append!(minima,Bounds_y[2])
                append!(maxima,Bounds_y[1])
            else
                append!(minima,Bounds_y[1])
                append!(maxima,Bounds_y[2])
            end
        end

        if Nvec[2] != 0
            Bounds_px = ([sqrt(2)*(-1 + 2/n *(j-1)),sqrt(2)*(-1 + 2/n *j)] - T[2]*[1.0,1.0])/Nvec[2]
            if Bounds_px[1] >= Bounds_px[2]
                append!(minima,Bounds_px[2])
                append!(maxima,Bounds_px[1])
            else
                append!(minima,Bounds_px[1])
                append!(maxima,Bounds_px[2])
            end
        end

        if Nvec[3] != 0
            Bounds_py = ([sqrt(2)*(-1 + 2/n *(i-1)),sqrt(2)*(-1 + 2/n *i)] - T[3]*[1.0,1.0])/Nvec[3]
            if Bounds_py[1] >= Bounds_py[2]
                append!(minima,Bounds_py[2])
                append!(maxima,Bounds_py[1])
            else
                append!(minima,Bounds_py[1])
                append!(maxima,Bounds_py[2])
            end
        end

        if abs(Nvec[1]) <= 1e-10 && abs(Nvec[2]) <= 1e-10 && abs(Nvec[3]) <= 1e-10
            println("uff")
            #println(List)
            amin = 0.0
            bmax = 0.0
        else
            amin = min(minima...)
            bmax = max(maxima...)
   
            if abs(amin) <= 1e-10 #find_zeros cant handle 1e-11 and lower
                amin = -1e-10
            end

            if abs(bmax) <= 1e-10
                bmax = 1e-10
            end

   
            Ham(t) = A*(x^2 + dot(vec(t),vec(t))) + B*(((vec(t)[2])^2 + (vec(t)[3])^2)*(2 - x^2 - dot(vec(t),vec(t))) + (x*vec(t)[3] - vec(t)[1]*vec(t)[2])^2)-E0 + C*vec(t)[3]*sqrt(Complex((2 - x^2 - dot(vec(t),vec(t)))))
            Ham2(t) = (A*(x^2 + dot(vec(t),vec(t))) + B*(((vec(t)[2])^2 + (vec(t)[3])^2)*(2 - x^2 - dot(vec(t),vec(t))) + (x*vec(t)[3] - vec(t)[1]*vec(t)[2])^2)-E0)^2 - C^2*vec(t)[3]^2*(2 - x^2 - dot(vec(t),vec(t)))
 
            point = find_zeros(Ham2,amin,bmax)
   
            ii = 1

            if Array{Float64,1}[] != point
        
                while Intersection === missing && ii<= length(point)
                    test = Ham(point[ii])
            
                    if imag(test) <1e-5 && abs(real(test)) < 1e-5 
                        TestIntersection = vec(real(test))
                        if sqrt(2)*(-1 + 2/n *(k-1)) <= TestIntersection[1] <= sqrt(2)*(-1 + 2/n *k) && sqrt(2)*(-1 + 2/n *(j-1)) <= TestIntersection[2] <= sqrt(2)*(-1 + 2/n *j) && sqrt(2)*(-1 + 2/n *(i-1)) <= TestIntersection[3] <= sqrt(2)*(-1 + 2/n *i)
               
                        Intersection = TestIntersection
                        end
                    end
                    ii+=1
                end
            end
        end
    end

    return Intersection #[y,px,py] nebo missing

end


# Intersection of line in box given by i,j,k coordinates with surface given by the vibron Hamiltonian
function Point_of_Int_line(X1::Array{Float64},X2::Array{Float64},par::Array{Float64})
    x,A,B,C,E0,i,j,k,n = par[1:9]

    Nvec=X2 - X1

    vec(t) = X1 + t*Nvec

    #urceni hranic t vazane na X,Y,Z
    minima = Float64[]
    maxima = Float64[]

    if abs(Nvec[1]) >= 1e-7
        Bounds_y = ([sqrt(2)*(-1 + 2/n *(k-1)),sqrt(2)*(-1 + 2/n *k)] - X1[1]*[1.0,1.0])/Nvec[1]

        if Bounds_y[1] >= Bounds_y[2]
            append!(minima,Bounds_y[2])
            append!(maxima,Bounds_y[1])
        else
            append!(minima,Bounds_y[1])
            append!(maxima,Bounds_y[2])
        end
    end

    if abs(Nvec[2]) >= 1e-7
        Bounds_px = ([sqrt(2)*(-1 + 2/n *(j-1)),sqrt(2)*(-1 + 2/n *j)] - X1[2]*[1.0,1.0])/Nvec[2]
        if Bounds_px[1] >= Bounds_px[2]
            append!(minima,Bounds_px[2])
            append!(maxima,Bounds_px[1])
        else
            append!(minima,Bounds_px[1])
            append!(maxima,Bounds_px[2])
        end
    end

    if abs(Nvec[3]) >= 1e-7
        Bounds_py = ([sqrt(2)*(-1 + 2/n *(i-1)),sqrt(2)*(-1 + 2/n *i)] - X1[3]*[1.0,1.0])/Nvec[3]
        if Bounds_py[1] >= Bounds_py[2]
            append!(minima,Bounds_py[2])
            append!(maxima,Bounds_py[1])
        else
            append!(minima,Bounds_py[1])
            append!(maxima,Bounds_py[2])
        end
    end
    koren = missing

    if abs(Nvec[1]) <= 1e-7 && abs(Nvec[2]) <= 1e-7 && abs(Nvec[3]) <= 1e-7
        println("uff")
    else
        amin = min(minima...) 
        bmax = max(maxima...)

        if abs(amin) <= 1e-10
        
            amin = 1e-10
        end
        if abs(bmax) <= 1e-10
        
            bmax = 1e-10
        end

        Ham(t) = A*(x^2 + dot(vec(t),vec(t))) + B*(((vec(t)[2])^2 + (vec(t)[3])^2)*(2 - x^2 - dot(vec(t),vec(t))) + (x*vec(t)[3] - vec(t)[1]*vec(t)[2])^2)-E0 + C*vec(t)[3]*sqrt(Complex((2 - x^2 - dot(vec(t),vec(t)))))
        Ham2(t) = (A*(x^2 + dot(vec(t),vec(t))) + B*(((vec(t)[2])^2 + (vec(t)[3])^2)*(2 - x^2 - dot(vec(t),vec(t))) + (x*vec(t)[3] - vec(t)[1]*vec(t)[2])^2)-E0)^2 - C^2*vec(t)[3]^2*(2 - x^2 - dot(vec(t),vec(t)))

        point = find_zeros(Ham2,amin,bmax)

    
        ii = 1

        if Array{Float64,1}[] != point
        
            while koren === missing && ii<= length(point)
                test = Ham(point[ii])
            
                if imag(test) <1e-5 && abs(real(test)) < 1e-5 
                    testkoren = vec(real(test))
                    if sqrt(2)*(-1 + 2/n *(k-1)) <= testkoren[1] <= sqrt(2)*(-1 + 2/n *k) && sqrt(2)*(-1 + 2/n *(j-1)) <= testkoren[2] <= sqrt(2)*(-1 + 2/n *j) && sqrt(2)*(-1 + 2/n *(i-1)) <= testkoren[3] <= sqrt(2)*(-1 + 2/n *i)
                    #println(koren)
                    koren = testkoren
                    end
                end
                ii+=1
            end
        end
    end
    return koren #[y,px,py] nebo missing

end

#
#   POINT PICKING
#
function Doubles(List1::Array{Array{Float64}},List2::Array{Array{Float64}},List3::Array{Array{Float64}},l1::Int64,l2::Int64,l3::Int64,par::Array{Float64})
    
    List1_Doubles =  Array{Array{Float64}}[]
    Point = missing                
    for ii in 1:l1
        for jj in 1:l1 - ii
            push!(List1_Doubles,[List1[ii],List1[ii + jj]])
        end
    end
    
    List_Supplements = vcat(List2,List3)
    Length_List1_Doubles = length(List1_Doubles)
    ii=1
    jj=1
    while  Point === missing && ii <= Length_List1_Doubles
        while  Point === missing && jj <= l2 + l3
                
            Point = Point_of_Intersection(vcat(List1_Doubles[ii],[List_Supplements[jj]]),par)
            jj+=1
        end
        ii+=1
    end
    
    return Point
end
    
    
function Triples(List1::Array{Array{Float64}},l1::Int64,par::Array{Float64})
    List1_Triples =  Array{Array{Float64}}[]
    Point = missing                
    for ii in 1:l1
        for jj in 1:l1 - ii
            for kk in 1:l1 - ii - jj
                push!(List1_Triples,[List1[ii],List1[ii + jj],List1[ii + jj + kk]])
            end
        end
    end
        
    Length_List1_Triples = length(List1_Triples)
    ii=1
    while  Point === missing && ii <= Length_List1_Triples
        Point = Point_of_Intersection(List1_Triples[1],par)
        ii+=1
    end
    
    return Point
end

function Points_and_Surface(n::Int64,par::Array{Float64},x::Float64; Plot::Array{Bool} = [false,false])

    E0,A,B,C = par[1:4]
    
    probb = 0 #Problematic points
    total = 0 #Total number of points
    S = 0.0

    if Plot[1] || Plot[2]
        pyplot()
        #pyplot(size = (2200,1000))
        PyPlot.pygui(true)
        fig = scatter()
    end

    plane1 = [InitialCondS_Px(x,sqrt(2)*(-1 + (i-1)*2/n),sqrt(2)*(-1 + (j-1)*2/n),E0,A,B,C) for i in 1:n+1, j in 1:n+1] #[y,py]   
      
    plane2 = [InitialCondS_Y(x,sqrt(2)*(-1 + (i-1)*2/n),sqrt(2)*(-1 + (j-1)*2/n),E0,A,B,C) for i in 1:n+1, j in 1:n+1] #[px,py]
        
    plane3 = [InitialCondS_Py(x,sqrt(2)*(-1 + (i-1)*2/n),sqrt(2)*(-1 + (j-1)*2/n),E0,A,B,C) for i in 1:n+1, j in 1:n+1] #[y,px]
    
    if Plot[2]

        for i in 1:n+1
            for j in 1:n+1
                if plane1[i,j] !== missing
                    
                ll = length(plane1[i,j])
                fig = scatter!( fill(sqrt(2)*(-1 + (i-1)*2/n),ll),plane1[i,j], fill(sqrt(2)*(-1 + (j-1)*2/n),ll), legend = false,c = "blue")
            
                end

                if plane2[i,j] !== missing
            
                    ll = length(plane2[i,j])
                    fig = scatter!(plane2[i,j], fill(sqrt(2)*(-1 + (i-1)*2/n),ll), fill(sqrt(2)*(-1 + (j-1)*2/n),ll), legend = false,c = "blue")
                    
                end

                if plane3[i,j] !== missing
            
                    ll = length(plane3[i,j])
                    fig = scatter!( fill(sqrt(2)*(-1 + (i-1)*2/n),ll),fill(sqrt(2)*(-1 + (j-1)*2/n),ll),plane3[i,j], legend = false,c = "blue")
                    
                end

            end
        end

    end

    points = Array{Array{Float64}}(undef, n, n, n)
    
    
    
    for i in 1:n #py
        for j in 1:n #px
            #println("making $i $j")
            for k in 1:n #y
                par = [x,A,B,C,E0,i,j,k,n]
    
                X = Array{Float64}[]
                Y = Array{Float64}[]
                Z = Array{Float64}[]
                lx = 0
                ly = 0
                lz = 0
                
                for ii in 1:2
                    for jj in 1:2
                        
                        if plane1[k + (ii-1),i + (jj-1)] !== missing #plane1[y,py]
                            
                            point = find(plane1[k + (ii-1),i + (jj-1)],sqrt(2)*(-1 + 2/n *(j-1)),sqrt(2)*(-1 + 2/n *j)) 
                            if point !== missing
                                
                                push!(X,[sqrt(2)*(-1 + 2/n *(k + (ii-1)-1)),point,sqrt(2)*(-1 + 2/n *(i + (jj-1)-1))])
                                lx += 1

                            end
                        end
                    end
                end
            
                for ii in 1:2
                    for jj in 1:2
                        if plane2[j + (jj-1),i + (ii-1)] !== missing #plane2[px,py]
                            point = find(plane2[j + (jj-1),i + (ii-1)],sqrt(2)*(-1 + (2/n)*(k-1)),sqrt(2)*(-1 + (2/n)*k)) 
                            if point !== missing

                                push!(Y,[point,sqrt(2)*(-1 + 2/n *(j + (jj-1)-1)),sqrt(2)*(-1 + (2/n)*(i + (ii-1) -1))])
                                ly += 1

                            end
                        end
                    end
                end
    
                for ii in 1:2
                    for jj in 1:2
                        if plane3[k + (ii-1),j + (jj-1)] !== missing #plane3[y,px]
                            point = find(plane3[k + (ii-1),j + (jj-1)],sqrt(2)*(-1 + 2/n *(i-1)),sqrt(2)*(-1 + 2/n *i)) 
                            if point !== missing
                                
                                push!(Z,[sqrt(2)*(-1 + 2/n *(k + (ii-1)-1)),sqrt(2)*(-1 + 2/n *(j + (jj-1)-1)),point])
                                lz += 1

                            end
                        end
                    end
                end
                
    
    
                
    
                    #roztrideni seznamu bodu podle delky
                if lx >= ly 
                    if lx >= lz
                        List1 = X
                        if ly >= lz
                            List2 = Y
                            List3 = Z
                        else 
                            List2 = Z
                            List3 = Y
                        end
                    else 
                        List1 = Z
                        List2 = X
                        List3 = Y
                    end
                elseif lz >= ly
                    List1 = Z
                    List2 = Y
                    List3 = X
                elseif lx >= lz
                    List1 = Y
                    List2 = X
                    List3 = Z
                else
                    List1 = Y
                    List2 = Z
                    List3 = X
                end
                Point  = missing
                l1 = length(List1)
                l2 = length(List2)
                l3 = length(List3)
    
                if lx + ly + lz >= 3
                    List = Array{Array{Float64}}[]
                    
                        
                    if l3 >= 1
                        for ii in 1:l1
                            for jj in 1:l2
                                for kk in 1:l3

                                    push!(List,[List1[ii],List2[jj],List3[kk]])
                                    
                                end 
                            end
                        end
    
                        Lenght_List = length(List)
                        iii = 1
                        while  Point === missing && iii <= Lenght_List
                            Point = Point_of_Intersection(List[iii],par)
                            iii+=1
                        end
                        
                    end
    
                    if Point === missing

                        #dvojice pro nejdelsi list        
                        Point = Doubles(List1,List2,List3,l1,l2,l3,par)
                        
                        if Point === missing
                            #dvojice pro druhy nejdelsi
                            Point = Doubles(List2,List1,List3,l2,l1,l3,par)
                            
                            if Point === missing

                                #dvojice pro treti nejdelsi
                                Point = Doubles(List3,List1,List2,l3,l1,l2,par)
                                
                                if Point === missing

                                    if l1 >=3
                                        Point = Triples(List1,l1,par)
                                    end

                                    if Point === missing

                                        if l2 >=3
                                            Point = Triples(List2,l2,par)
                                        end
    
                                        if Point === missing
    
                                            if l3 >=3
                                                Point = Triples(List3,l3,par)
                                            end
                                            
                                        end
                                    end
                                end
    
                            end
                        end
                    end 
    
                    if Point === missing
                        
                        probb+=1
                        points[k,j,i] = append!([x],List1[1])
                        S_kji = Surface_linear(List1[1],par)
                        append!(points[k,j,i],S_kji)
                        S += S_kji
                        total +=1
                    else 
                
                        points[k,j,i] = append!([x],Point)
                        S_kji = Surface_linear(Point,par)
                        append!(points[k,j,i],S_kji)
                        S += S_kji
                        total +=1
                    end
    
                elseif l1 + l2 + l3 == 1
                    points[k,j,i] = append!([x],List1[1])
                    S_kji = Surface_linear(List1[1],par)
                    append!(points[k,j,i],S_kji)
                    S += S_kji
                    total +=1
    
                elseif l1 + l2 == 2

                    if l1 == 2
                        Point = Point_of_Int_line(List1[1],List1[2],par)
                    else
                        Point = Point_of_Int_line(List1[1],List2[1],par)
                    end
                        
                    if Point !== missing 
                        points[k,j,i] = append!([x],Point)
                        S_kji = Surface_linear(List1[1],par)
                        append!(points[k,j,i],S_kji)
                        S += S_kji
                        total +=1
                    else
                        points[k,j,i] = [NaN]
                    end
    
              
                else
                    points[k,j,i] = [NaN]
                end
    
                if Plot[1] && !isnan(points[k,j,i][1])
                    fig = scatter!( [points[k,j,i][2]],[points[k,j,i][3]],[points[k,j,i][4]], legend = false,c = "green")
                    #println(E_vibron(points[k,j,i][1:4],[A,B,C]))
                end
    
            end
        end
    end
  
    if Plot[1]||Plot[2]
        display(fig)
    end
    return points,S,total

end
#println(Points_and_Surface(50,[0.0,0.1,-0.8,0.25],0.0,Plot = [false,true]))


#
#
#   EVALUATION OF LYAPUNOV EXPONENT
#
#

function tr_condition(t,u,integrator)
    true
end


function affect_norm!(integrator)
    
    a = integrator.p[4][1]
    integrator.p[4][1] =integrator.t[end]
    len = length(integrator.p[3])
    
    n = norm(integrator.u[5:8])
    
    if len == 0
        push!(integrator.p[3],log(n)/integrator.t[end])

    
    else
        l = (integrator.p[3][len]*a + log(n))/integrator.t[end]
        
        if len == 400
            
            integrator.p[3][1:399] = integrator.p[3][2:400]
            integrator.p[3][400] = l
           
             
            Var = var(integrator.p[3])
            Mean = mean(integrator.p[3])
            
            if abs(Var/Mean) < 1e-5 || abs(Mean) < 1e-5
                #return terminate!(integrator) 
            end
    
        else
                    
            push!(integrator.p[3],l)
        end
    end
    integrator.u[5:8] /= n
        
end

function Vibron!(du,u,p,t)
   
    A, B, C = p[1][2:4]
    x,y,px,py = u[1:4]


    Σ = x^2 + y^2 + px^2 + py^2

    if Σ > 2
        println("uff")
        println(x,y,px,py)
              
        println(Σ)
        @error("s negative!")
        
    else

        M = x*px + y*py
        D = px^2 + py^2 -1
        

        dCn = 2*[px, py, -x, -y]
          
        dCW = 2*[-(2*px*D+x*M),-(2*py*D+y*M),px*M,py*M]
        dCLx = [-px*py,-(-(2-Σ)+py^2),py*x,py*y]/((2 - Σ)^(1/2))
      
        du[1:4] = A*dCn + B*dCW + C*dCLx
   
    HCn = 2*[ 0 0 1 0;
           0 0 0 1;
          -1 0 0 0;
          0 -1 0 0]

    HCW = [ -2*(2*px*x + py*y)          -2*py*x           -2*(-2 + 6*px^2 + 2*py^2 + x^2)      -2(4*px*py + x*y);
                -2*px*y           -2*(2*py*y + px*x)       -2(4*px*py + x*y)           -2*(-2 + 6*py^2 + 2*px^2 + y^2);
                  2*px^2                 2*px*py             2*(2*px*x + py*y)                       2*px*y;
                 2*px*py                  2*py^2                  2*py*x                        2*(px*x + 2*py*y)]   

    HCLx =     [ -px*py*x                     -px*py*y                py*(-2 +Σ -px^2)       px*(-2 +Σ -py^2);
           x*(-2 +Σ -py^2)          y*(-2 +Σ -py^2)                px*(-2 +Σ -py^2)       py*(-6 +3*Σ -py^2);
           -py*(-2 +Σ -x^2)                py*x*y                       px*py*x           -x*(-2 +Σ -py^2);
               py*x*y                  -py*(-2 +Σ -y^2)                 px*py*y            -y*(-2 +Σ -py^2)]/(2 - Σ)^(3/2)

    du[5:8] = (A*HCn + B*HCW + C*HCLx)*u[5:8]

    end

end


function g(resid, u, p, t)
    resid[1] = E_vibron(u,p[1][2:4]) - p[1][1]
    resid[2] = 0

end


function intersection_x(u,t,integrator)
    #u[Int64(integrator.p[2])] - integrator.p[5]
    u[1] - integrator.p[5][1]
 
end

function affect!(integrator)
     
end



function CheckDomain(u,p,t)
   x, y, px, py = u[1:4]
   
   s = 2.0 - (x^2 + y^2 + px^2 + py^2)
   return s <= 0
end



function Lyapunov!(Point::Array{Float64},par::Array{Float64},x::Float64)
    E0,A,B,C = par[1:4]
 
    u0 = Point
    #println(Point)
    #println(E_vibron(Point,par[2:4]))
    v = [rand(),rand(),rand(),rand()]
    v = v/norm(v)
    append!(u0,v)

    p = ([E0,A,B,C],[1],Float64[],[0.0],[x])
    tspan=(0.0,1e5)
        
    prob = ODEProblem(Vibron!, u0, tspan, p)   

    
    cb_intersect = ContinuousCallback(intersection_x, affect!;save_positions = (true,false))
    cb_norm = DiscreteCallback(tr_condition, affect_norm!;save_positions = (false,false))
    cb_man = ManifoldProjection(g; save = false)

    cb = CallbackSet(cb_intersect,cb_norm)#,cb_man)
    #println("tu1")    
    time = @elapsed solution = try
                                 solve(prob,Vern9(), callback = cb, abstol = 1e-10,reltol = 1e-10, maxiters = 1e7, save_everystep=false, save_start=true, save_end=false, isoutofdomain = CheckDomain)       
                            catch
                                 missing
                            end
#println("tu")
    #time = @elapsed solution = solve(prob,Vern9(), callback = cb, abstol = 1e-10,reltol = 1e-10, maxiters = 1e7, save_everystep=false, save_start=true, save_end=false, isoutofdomain = CheckDomain)       
    #println(solution)
    #println(solution.retcode)
    if solution !== missing# && solution.retcode != :DtNaN 
        if solution.retcode != :DtNaN 
            #println("in")                   
            lyapunov = mean(p[3])
            lv = var(p[3])
    
        #println("Time = $time, Lyapunov = $lyapunov ± $lv")
        #println(solution.retcode)

            Ax = [a[1] for a in solution.u]
            Ay =[a[2] for a in solution.u]
            Apx = [a[3] for a in solution.u]
            Apy = [a[4] for a in solution.u]
            if !isnan(lyapunov)
                return solution.u, lyapunov#[[Ax[ii],Ay[ii],Apx[ii],Apy[ii]] for ii in 1:length(Ax)],lyapunov
            else 
                return missing, missing
            end
        else
            return missing, missing
        end

    else
        #println("Time = $time")
        #println(solution.retcode)
        return missing, missing
    end

end