#A set of 2D interacting particles. The interaction is modelled using a Harmonic Potential.
#The particles are initalized from auniform grid. 
using Plots
plotlyjs()
##
function V_H(xl,yl)# Defing the interaction between particles to be harmonic 
    x = xl[1] - xl[2]
    y = yl[1] - yl[2]
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r < r_cut
    	v = 0.5*K*(r-r_eq)^2#harmonic potential in 2 dimensions
    else
    	v = 0.5*K*(r_cut-r_eq)^2
    end
    return v
end
##
function V_H_plot(xl,yl)# Defing the interaction between particles to be harmonic 
    x = xl
    y = yl
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r < r_cut
    	v = 0.5*K*(r-r_eq)^2#harmonic potential in 2 dimensions
    else
    	v = 0.5*K*(r_cut-r_eq)^2
    end
    return v
end
##
#Harmonic intearction is assumed between particles
function Fint_H(xl,yl,k=K)
    x = xl[1] - xl[2]
    y = yl[1] - yl[2]
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r < r_cut
    	Fx = -k*x*(1-r_eq/r)
    	Fy = -k*y*(1-r_eq/r)
    else
    	Fx, Fy = 0, 0
    end
    return Fx,Fy
end
##
function Fx_H(xl,yl,k=K)
    x = xl
    y = yl
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r < r_cut
    	Fx = -k*x*(1-r_eq/r)
    	Fy = -k*y*(1-r_eq/r)
    else
    	Fx, Fy = 0, 0
    end
    return Fx
end
function Fy_H(xl,yl,k=K)
    x = xl
    y = yl
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r < r_cut
    	Fx = -k*x*(1-r_eq/r)
    	Fy = -k*y*(1-r_eq/r)
    else
    	Fx, Fy = 0, 0
    end
    return Fy
end
##
function Run(Xinit=0,Yinit=0,Vxinit=0,Vyinit=0;T=1:0.1:10,np=10)
    #INPUT: initial cordinates and velocities and the time axis and number of particles
    #RETURNS: The X,Y as a function of time. Total energy as a function of time
    lent = length(T)
    X = zeros(lent,np)
    Y = zeros(lent,np)
    E = zeros(lent)
    Temp = zeros(lent)#the temperature i.e, the average kinetic energy
    PE = zeros(lent)
    xnew = Xinit
    ynew = Yinit
    vxnew = Vxinit
    vynew = Vyinit
    for i in 1:lent
        #fxy = F.(xnew,[Xinit[3],Xinit[6]]ynew)
        #fxt = [ f[1] for f in fxy]
        #fyt = [ f[2] for f in fxy]
          #parameter for boundary
        xinit = copy(xnew)
        #
        fxt = zeros(np)
        fyt = zeros(np)
        for i in 1:length(xnew)-1
            for j in i+1:length(ynew)
                fxint,fyint = Fint_H([xnew[i],xnew[j]],[ynew[i],ynew[j]])
                fxt[i] += fxint
                fyt[i] += fyint
                fxt[j] += -fxint
                fyt[j] += -fyint
            end
        end
        xnew = xnew .+ vxnew*dt .+ 0.5*fxt*dt*dt./m
        ynew = ynew .+ vynew*dt .+ 0.5*fyt*dt*dt./m
        #fxyn = F.(xnew, ynew)
        #fxtn = [ f[1] for f in fxyn]
        #fytn = [ f[2] for f in fxyn]
        fxtn = zeros(np)
        fytn = zeros(np)
        PE_elec = 0 #Potential energy due to electronic interactions
        for i in 1:length(xnew)-1
            for j in i+1:length(ynew)
                fxint,fyint = Fint_H([xnew[i],xnew[j]],[ynew[i],ynew[j]])
                fxtn[i] += fxint
                fytn[i] += fyint
                fxtn[j] += -fxint
                fytn[j] += -fyint
                PE_elec += V_H([xnew[i],xnew[j]],[ynew[i],ynew[j]])
            end
        end
        vxnew = vxnew .+ 0.5*(fxt .+ fxtn)./m*dt
        vynew = vynew .+ 0.5*(fyt .+ fytn)./m*dt
        X[i,:] = xnew
        Y[i,:] = ynew
        PE[i] = PE_elec
        KE = sum( ( 0.5*m.*(vxnew.^2 + vynew.^2)) )
        E[i] =  KE + PE_elec
        Temp[i] = KE/np
        for i in 1:np
            #global box in x
            if (abs(xnew[i]) > bl)
                vxnew[i] = - vxnew[i]
            #left subbox boundary
            elseif (((xinit[i]+tl)*(xnew[i]+tl))<0 && abs(ynew[i])>tw)
                vxnew[i] = - vxnew[i]
                xnew[i] = xinit[i]
            elseif (((xinit[i]-tl)*(xnew[i]-tl))<0 && abs(ynew[i])>tw)
                vxnew[i] = - vxnew[i]
                xnew[i] = xinit[i]
            end
            #global box
            if (abs(xnew[i]) > tl)
                if (abs(ynew[i]) > bl)
                    vynew[i] = - vynew[i]
                end
            else
                if (abs(ynew[i]) > tw)
                    vynew[i] = - vynew[i]
                end
            end
        end
    end
    return X,Y,E,Temp,PE
end
##Interaction Parameters
np = 20
K = 3
r_eq = 2
r_cut = 3#2.8 #the cut-off distance for interactions
##box parameters
bl = 20 # half-box length
tl = 5  # tunnel lengthplot(T,Te)
tw = 5  # tunnel width
#animation parameter
frames = 500
##
print("Setting Initial Conditions")
dt =0.0005
T = 0:dt:30
lt = length(T)
m = ones(np)
##  Grid Initialization of position coordinates 
set_x =  -19:3:-9#LinRange(-19.9,-5.1,4)
set_y =  6:3:19#LinRange(-19.9,19.9,5)
Xinit = []#zeros(np)
Yinit = []#zeros(np)
#particle index
for x in set_x
    for y in set_y
        append!(Xinit,x)
        append!(Yinit,y)
    end
end
##
Vxinit = 8*(rand(np).-0.5)
Vyinit = 5*(rand(np).-0.5)
##
print("\nRunning Simulation") 
Xt,Yt,Et,Te,PE=Run(Xinit,Yinit,Vxinit,Vyinit;T=T,np=np)
print("\nEnergy Error:",maximum(Et)-minimum(Et))
##
print("\n Plotting")
plot(T,Et)
savefig("Total_Energy.png")
plot(T,Te)
savefig("Temperature.png")
##
#plot(T,PE)
##
surface(-r_cut-0.5:0.1:r_cut+0.5,-r_cut-0.5:0.1:r_cut+0.5,V_H_plot)
savefig("V.png")
##
N_step_frame =  Int(ceil(lt/frames))#No of steps per frame
print("\nAnimating")
plot()
anim = @animate for j in 1:frames
    i = Int(ceil(j*lt/frames))
    if j > 4
        plot([-tl,-bl,-bl,-tl],[bl,bl,-bl,-bl], color = :black,linewidth=3)
        plot!([tl,bl,bl,tl],[bl,bl,-bl,-bl], color = :black,linewidth=3)
        plot!([-tl,-tl],[-bl,-tw], color = :black,linewidth=3,ticks = false)
        plot!([-tl,-tl],[bl,tw], color = :black,linewidth=3)
        plot!([tl,tl],[-bl,-tw], color = :black,linewidth=3,ticks = false)
        plot!([tl,tl],[bl,tw], color = :black,linewidth=3)
        plot!([-tl,tl],[-tw,-tw], color = :black,linewidth=3)
        plot!([-tl,tl],[tw,tw], color = :black,linewidth=3)
        scatter!([Xt[i,:]],[Yt[i,:]],xlims=(-25,25),ylims=(-25,25),legend=false,color=:blue,markersize=2)
        plot!(Xt[(i-4*N_step_frame):i,:],Yt[(i-4*N_step_frame):i,:],legend=false,color=:blue,linewidth=1)
    else
        plot([-tl,-bl,-bl,-tl],[bl,bl,-bl,-bl], color = :black,linewidth=3)
        plot!([tl,bl,bl,tl],[bl,bl,-bl,-bl], color = :black,linewidth=3)
        plot!([-tl,-tl],[-bl,-tw], color = :black,linewidth=3,ticks = false)
        plot!([-tl,-tl],[bl,tw], color = :black,linewidth=3)
        plot!([tl,tl],[-bl,-tw], color = :black,linewidth=3,ticks = false)
        plot!([tl,tl],[bl,tw], color = :black,linewidth=3)
        plot!([-tl,tl],[-tw,-tw], color = :black,linewidth=3)
        plot!([-tl,tl],[tw,tw], color = :black,linewidth=3)
        scatter!([Xt[i,:]],[Yt[i,:]],xlims=(-25,25),ylims=(-25,25),legend=false,color=:blue,markersize=2)
        plot!(Xt[1:i,:],Yt[1:i,:],legend=false,color=:blue,linewidth=1)
    end
    #plot!([Xt[i:i,:]],[Yt[i-1000:i,:]],legend=false,dpi=500,color=:blue)
end
gif(anim,"Tunnel_box.gif",fps=30)
##
print("\n Done")
##
#surface(-6:0.01:6,-6:0.01:6,Fy_H)
