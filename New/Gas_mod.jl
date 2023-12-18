#A set of 2D interacting particles. The interaction is modelled using a Harmonic Potential.
#The particles are initalized radomly. 
##
using Plots
plotlyjs()
##
function V_H(xl,yl,k=K)# Defing the interaction between particles to be harmonic 
    x = xl[1] - xl[2]
    y = yl[1] - yl[2]
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r <= r_cut
    	v = 0.5*k*(r-r_eq)^2#harmonic potential in 2 dimensions
    else
    	v = 0.5*k*(r_cut-r_eq)^2
    end
    return v
end
##
#Harmonic intearction is assumed between particles
function Fint_H(xl,yl,k=K)
    x = xl[1] - xl[2]
    y = yl[1] - yl[2]
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    if r <= r_cut
    	Fx = -k*x*(1-r_eq/r)
    	Fy = -k*y*(1-r_eq/r)
    else
    	Fx, Fy = 0, 0
    end
    
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
        #fxy = F.(xnew,ynew)
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
        E[i] = sum( ( 0.5*m.*(vxnew.^2 + vynew.^2)) ) + PE_elec
        Temp[i] = sum(0.5*m.*(vxnew.^2 + vynew.^2))/np
        ##Reflections and boundary
        for i in 1:np
            #global box in x
            if (abs(xnew[i]) > 20)
                vxnew[i] = - vxnew[i]
            #left subbox boundary
            elseif (((xinit[i]+tl)*(xnew[i]+tl))<0 && abs(ynew[i])>wl)
                vxnew[i] = - vxnew[i]
                xnew[i] = xinit[i]
            elseif (((xinit[i]-tl)*(xnew[i]-tl))<0 && abs(ynew[i])>wl)
                vxnew[i] = - vxnew[i]
                xnew[i] = xinit[i]
            end
            #global box
            if (abs(xnew[i]) > tl)
                if (abs(ynew[i]) > 20)
                    vynew[i] = - vynew[i]
                end
            else
                if (abs(ynew[i]) > wl)
                    vynew[i] = - vynew[i]
                end
            end
        end
    end
    return X,Y,E,Temp
end
##Parameters
np = 20
K = 1
r_eq = 5
r_cut = 7#the cut-pff distance for interactions
wl = 5#3
tl = 5
#animation parameter
frames = 500
##
print("Setting Initial Conditions")
dt =0.0003
T = 0:dt:30
lt = length(T)
m = ones(np)
Xinit = 14*(rand(np).-0.5).-12.5
Yinit = 39*(rand(np).-0.5)
Vxinit = 20*(rand(np).-0.5)
Vyinit = 20*(rand(np).-0.5)

##
print("\nRunning Simulation")
Xt,Yt,Et,Te=Run(Xinit,Yinit,Vxinit,Vyinit;T=T,np=np)
print("\nEnergy Error:",maximum(Et)-minimum(Et))
##
print("\n Plotting")
plot(T,Et)
savefig("Total_Energy.png")
plot(T,Te)
savefig("Temperature.png")
##
N_step_frame =  Int(ceil(lt/frames))#No of steps per frame
print("\nAnimating")
plot()
anim = @animate for j in 5:frames
    i = Int(ceil(j*lt/frames))
    plot([-5,-20,-20,-5],[20,20,-20,-20], color = :black,linewidth=3)
    plot!([5,20,20,5],[20,20,-20,-20], color = :black,linewidth=3)
    plot!([-5,-5],[-20,-wl], color = :black,linewidth=3,ticks = false)
    plot!([-5,-5],[20,wl], color = :black,linewidth=3)
    plot!([5,5],[-20,-wl], color = :black,linewidth=3,ticks = false)
    plot!([5,5],[20,wl], color = :black,linewidth=3)
    plot!([-5.1,5.1],[-wl,-wl], color = :black,linewidth=3)
    plot!([-5.1,5.1],[wl,wl], color = :black,linewidth=3)
    scatter!([Xt[i,:]],[Yt[i,:]],xlims=(-25,25),ylims=(-25,25),legend=false,color=:blue,markersize=2)
    plot!(Xt[(i-4*N_step_frame):i,:],Yt[(i-4*N_step_frame):i,:],legend=false,color=:blue,linewidth=1)
    #plot!([Xt[i:i,:]],[Yt[i-1000:i,:]],legend=false,dpi=500,color=:blue)
end
##
gif(anim,"Tunnel_box.gif",fps=30)
##
print("\n Done")
