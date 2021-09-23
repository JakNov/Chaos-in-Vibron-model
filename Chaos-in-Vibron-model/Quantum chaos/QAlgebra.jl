using QuantumAlgebra

using Plots
using PyCall

@boson_ops tp
@boson_ops tm

@boson_ops s


#
#   zkouška balíčku QuantumAlgebra - přímý pokus o algebraický výpočet 
#   časového vývoje
#

#   zavedení stavu pomocí kreačních a anihilačních operátorů
function State(p::Int64,m::Int64,N::Int64)

    State  = tpdag()^p * tmdag()^m * sdag()^(N-m-p)
    return State
end


ξ = 0.5

n = normal_form(tpdag()*tp() + tmdag()*tm()) 
l = normal_form(tpdag()*tp() - tmdag()*tm())

Dp = normal_form(sqrt(2)*(tpdag()*s() - sdag()*tm()))
Dm = normal_form(sqrt(2)*(-tmdag()*s() + sdag()*tp()))


Rp = normal_form(sqrt(2)*(tpdag()*s() + sdag()*tm()))
Rm = normal_form(sqrt(2)*(tmdag()*s() + sdag()*tp()))

H = normal_form(n + 1/2*(Dp*Dm + Dm*Dp) +l*l)

state = State(1,0,3)

otoc = [1.0]

O = normal_form(comm(n,H))
val  = julia_expression(vacExpVal(O,state))

if abs(val) < 0.5
    val = 0.0
end
    
append!(otoc,val)


# vypocet koeficientu v sume
for i in 2:10
    global O = normal_form(comm(O,H))
    #println(O)
    val  = julia_expression(vacExpVal(O,state))

    if abs(val) < 0.5 
        val = 0.0
    end 
        
    append!(otoc,val/factorial(i))
    println(i)
end

#vypocet evol. operatoru pomoci sumy - rychly rust - moc kolem desateho kroku
function f(t,coeff)
    f = 0
    for i in 1:length(coeff)
        f = f + im^(i-1)*coeff[i] * t^(i-1)
    end
    return f
end 

println(f(1.0, otoc))

pyplot()
PyPlot.pygui(true)

time = [i/100 for i in 1:1000]
plot(time,  [real(f(time[i],otoc)*f(time[i],otoc)') for i in 1:1000])
#normal_form(comm(comm(comm(n,H),H),H))