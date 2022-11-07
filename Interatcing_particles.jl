#using Makie
using Plots
plotlyjs()
##
V(x,y) = 1*exp((x^2 + y^2)^0.5 - 20)#The global potential
function F(x,y)#The Glocal Force in which the gas evolves
    Fx = -x*exp((x^2 + y^2)^0.5 - 20)*(x^2+y^2)^-0.5
    Fy = -y*exp((x^2 + y^2)^0.5 - 20)*(x^2+y^2)^-0.5
    return Fx,Fy
end
function Fint(xl,yl)
    x1, y1 = xl[1], yl[1]
    x2, y2 = xl[2], yl[2]
    wx,wy=5,5
    kx=500
    ky=500
    Fx = kx*exp(-wx*(x2 - x1)^2)*sign(x1-x2)
    Fy = ky*exp(-wy*(y2 - y1)^2)*sign(y1-y2)
    return Fx, Fy
end
function V_elec(xl,yl)#electrostatic potential
    x1, y1 = xl[1], yl[1]
    x2, y2 = xl[2], yl[2]
    ϕ = 10#Some constant of repulsion/attraction
    ρ = ((x1-x2)^2 + (y1-y2)^2)#Distance between the particles squared
    return ϕ/ρ^0.5
end
#An electrostatic interaction is considered for this new system
function Fint_ele(xl,yl)
    x1, y1 = xl[1], yl[1]
    x2, y2 = xl[2], yl[2]
    ϕ = 10#Some constant of repulsion/attraction
    ρ = ((x1-x2)^2 + (y1-y2)^2)#Distance between the particles squared
    Fx = ϕ*(x1-x2)/ρ^1.5
    Fy = ϕ*(y1-y2)/ρ^1.5
    return Fx,Fy
end
##
function Run(Xinit=0,Yinit=0,Vxinit=0,Vyinit=0;T=1:0.1:10,np=10)
    #INPUT: initial cordinates and velocities and the time axis and number of particles
    #RETURNS: The X,Y as a function of time. Total energy as a function of time
    X = zeros(length(T),np)
    Y = zeros(length(T),np)
    E = zeros(length(T))
    Temp = zeros(length(T))#the temperature i.e, the average kinetic energy
    xnew = Xinit
    ynew = Yinit
    vxnew = Vxinit
    vynew = Vyinit
    for i in 1:length(T)
        fxy = F.(xnew,ynew)
        fxt = [ f[1] for f in fxy]
        fyt = [ f[2] for f in fxy]
        for i in 1:length(xnew)-1
            for j in i+1:length(ynew)
                fxint,fyint = Fint_ele([xnew[i],xnew[j]],[ynew[i],ynew[j]])
                fxt[i] += fxint
                fyt[i] += fyint
                fxt[j] += -fxint
                fyt[j] += -fyint
            end
        end
        xnew = xnew .+ vxnew*dt .+ 0.5*fxt*dt*dt./m
        ynew = ynew .+ vynew*dt .+ 0.5*fyt*dt*dt./m
        fxyn = F.(xnew, ynew)
        fxtn = [ f[1] for f in fxyn]
        fytn = [ f[2] for f in fxyn]
        PE_elec = 0 #Potential energy due to electronic interactions
        for i in 1:length(xnew)-1
            for j in i+1:length(ynew)
                fxint,fyint = Fint_ele([xnew[i],xnew[j]],[ynew[i],ynew[j]])
                fxtn[i] += fxint
                fytn[i] += fyint
                fxtn[j] += -fxint
                fytn[j] += -fyint
                PE_elec += V_elec([xnew[i],xnew[j]],[ynew[i],ynew[j]])
            end
        end
        vxnew = vxnew .+ 0.5*(fxt .+ fxtn)./m*dt
        vynew = vynew .+ 0.5*(fyt .+ fytn)./m*dt
        X[i,:] = xnew
        Y[i,:] = ynew
        E[i] = sum( (V.(xnew,ynew) + 0.5*m.*(vxnew.^2 + vynew.^2)) ) + PE_elec
        Temp[i] = sum(0.5*m.*(vxnew.^2 + vynew.^2))/np
    end
    return X,Y,E,Temp
end
##
dt =0.002
T = 0:dt:40
lt = length(T)
np = 50
m = ones(np)
Xinit = 20*(rand(np).-0.5)
Yinit = 20*(rand(np).-0.5)
Vxinit = 10*(rand(np).-0.5)
Vyinit = 10*(rand(np).-0.5)
##
Xt,Yt,Et,Te=Run(Xinit,Yinit,Vxinit,Vyinit;T=T,np=np)
##
plot(Et)
savefig("Total_Energy.png")
print(maximum(Et)-minimum(Et))
##
plot(T,Te)
savefig("Temperature.png")
##
tjn=1:50
plot(Xt[:,tjn],Yt[:,tjn])
##
#plot([-5,-5],[-5,5],linecolor=["black"])
#plot!([5,5],[-5,5],linecolor=["black"])
#plot!([-5,5],[5,5],linecolor=["black"])
#plot!([-5,5],[-5,-5],linecolor=["black"])
plot()
anim = @animate for j in 1:400
    i = Int(ceil(j*lt/400))
    scatter([Xt[i,:]],[Yt[i,:]],xlims=(-50,50),ylims=(-50,50),legend=false,color=:blue)
    #plot!([X[i-1000:i,:]],[Y[i-1000:i,:]],legend=false,dpi=500,color=:blue)
end
##
gif(anim,"Repulsieve_electron_gas_50.gif",fps=20)
