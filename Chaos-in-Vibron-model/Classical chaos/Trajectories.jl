

include("Functions.jl")

#
#   evolution of two close trajectories; returns their distance in time
#

function CloseTrajectories!(Point::Array{Float64},par::Array{Float64})

    x,y,py = Point
    E0,A,B,C = par[1:4]
    px = InitialCondS_Px(x,y,py,E0,A,B,C)
    Point = [x,y,px[1],py]
    u0 = [x,y,px[1],py]

    println([x,y,px[1],py])
    pyplot()
    PyPlot.pygui(true)

    v = [rand(),rand(),rand(),rand()]
    v = v/norm(v)
    append!(u0,v)

    p = ([E0,A,B,C],[1],Float64[],[0.0],[x])
    tspan=(0.0,1e4)
            
    prob = ODEProblem(Vibron!, u0, tspan, p)   

        
    cb_intersect = ContinuousCallback(intersection_x, affect!;save_positions = (true,false))
    cb_norm = DiscreteCallback(tr_condition, affect_norm!;save_positions = (false,false))
    cb_man = ManifoldProjection(g; save = false)

    cb = CallbackSet(cb_intersect,cb_norm)#,cb_man)
        #println("tu1")    
    time = @elapsed solution = try
                                    solve(prob,Vern9(), callback = cb, abstol = 1e-10,reltol = 1e-10, maxiters = 1e7, save_everystep=true, save_start=true, save_end=true, isoutofdomain = CheckDomain)       
                                catch
                                    missing
                                end
    println("Time = $time")
    println(solution.retcode)
                                
    plot(solution,vars=(1,3),denseplot=true)
    lyapunov1 = mean(p[3])
    name = "$Point lyp = $lyapunov1"
                                
    println(name)
    png(name)

    name = "Lyap1 $Point lyp = $lyapunov1"
    plot(p[3])                
    println(name)
    png(name)
    Sol1 = solution.u
    Lyap1 = p[3]
    
    x = x +1e-6
    px = InitialCondS_Px(x,y,py,E0,A,B,C)
    Point = [x,y,px[1],py]

    u0 = [x,y,px[1],py]

    println([x,y,px[1],py])
    pyplot()
    PyPlot.pygui(true)

    v = [rand(),rand(),rand(),rand()]
    v = v/norm(v)
    append!(u0,v)

    p = ([E0,A,B,C],[1],Float64[],[0.0],[x])
    tspan=(0.0,1e4)
            
    prob = ODEProblem(Vibron!, u0, tspan, p)   

        
    cb_intersect = ContinuousCallback(intersection_x, affect!;save_positions = (true,false))
    cb_norm = DiscreteCallback(tr_condition, affect_norm!;save_positions = (false,false))
    cb_man = ManifoldProjection(g; save = false)

    cb = CallbackSet(cb_intersect,cb_norm)#,cb_man)
        #println("tu1")    
    time = @elapsed solution = try
                                    solve(prob,Vern9(), callback = cb, abstol = 1e-10,reltol = 1e-10, maxiters = 1e7, save_everystep=true, save_start=true, save_end=true, isoutofdomain = CheckDomain)       
                                catch
                                    missing
                                end
    println("Time = $time")
    println(solution.retcode)
                                
    plot(solution,vars=(1,3),denseplot=true)
    lyapunov2 = mean(p[3])
    name = "$Point lyp = $lyapunov2"
                                
    println(name)
    png(name)

    name = "Lyap2 $Point lyp = $lyapunov2"
    plot(p[3])                
    println(name)
    png(name)
    Sol2 = solution.u
    Lyap2 = p[3]

    Distance = []
    ln = min(length(Sol1[:]),length(Sol2[:]))
    for i in 1:ln
        d = sum((Sol1[i] - Sol2[i]).^2)
        #println(d)
        append!(Distance,[d])
    end

    name = "Distance $Point $par $lyapunov1 $lyapunov2"
    plot(Distance)           
    println(name)
    png(name)
end              

#CloseTrajectories!([0.1,-0.2,0.2],[0.0,0.1,-0.8,0.0])
#CloseTrajectories!([0.1,-0.2,0.2],[0.0,0.1,-0.8,0.1])
#CloseTrajectories!([0.1,-0.2,0.2],[0.0,0.1,-0.8,0.2])
#CloseTrajectories!([0.1,-0.2,0.2],[0.0,0.1,-0.8,0.3])
#CloseTrajectories!([0.1,-0.2,0.2],[0.0,0.1,-0.8,0.4])
#CloseTrajectories!([0.1,-0.2,0.2],[0.0,0.1,-0.8,0.5])

