
# Fits a 2D Keplerian to a set of data
using PyPlot
using Optim

# Kepler solver:
include("../../src/kepler.jl")

function return_xy(t,elements)
  P,ecc,tp,semi,omega = elements
  if ecc >= 1.0
    return 0.0,0.0
  end
  # First, compute M:
  M = (t-tp)*2pi/P
  # Now, compute eccentric anomaly:
  ekep = ekepler(M,ecc)
  # Now compute xp & yp (the planar values aligned with pericenter):
  xp = semi*(cos(ekep)-ecc)
  yp = semi*sqrt(1-ecc^2)*sin(ekep)
  # Now rotate:
  x = xp*cos(omega)-yp*sin(omega)
  y = xp*sin(omega)+yp*cos(omega)
return x,y
end

function fit_keplerian(t,x,y,elements0)

  # elements are given by: [P,ecc,tp,semi,omega]
  # Chi-square function:
  function chi_square(elements)
    # Call model:
    chisq = 0.0
    for i=1:length(t)
      xmod,ymod = return_xy(t[i],elements)
      chisq += (x[i]-xmod)^2+(y[i]-ymod)^2
    end
#    println("elements: ",elements," chi^2: ",chisq)
  return chisq
  end
   
  # Optimize the fit:
  result =  optimize(chi_square,elements0,LBFGS(),autodiff=:forward)
#  result =  optimize(chi_square,elements0,LBFGS())
  elements_best = result.minimizer
  xmod = zeros(t)
  ymod = zeros(t)
  for i=1:length(t)
     xmod[i],ymod[i] = return_xy(t[i],elements_best)
  end
return xmod,ymod, result
end

data_hr  = readdlm("test_output_highres.txt")
elements_all = readdlm("elements.txt",',')
semi = median(sqrt.(data_hr[:,5].^2 .+data_hr[:,7].^2))
period = elements_all[2,2]
elements_best = zeros(5)
chi_best = Inf
chisq = zeros(20,20)
tobs = data_hr[:,1]; xobs  = data_hr[:,5]-data_hr[:,2]; yobs = -(data_hr[:,7]-data_hr[:,4])
xmod = zeros(tobs)
ymod = zeros(tobs)
k=1
for j=1:20
  tp = elements_all[2,3]+(j-0.5)/20*elements_all[2,2]
  M = (tobs[1]-tp)*2*pi/period
  omega = atan2(-(yobs[1]*cos(M)-xobs[1]*sin(M))/semi,(yobs[1]*sin(M)+xobs[1]*cos(M))/semi)
  for k=1:20
    omega = (k-0.5)/20*2*pi
    elements = [period,sqrt(elements_all[2,4]^2+elements_all[2,5]^2),tp,semi,omega]
  # P,ecc,tp,semi,omega = elements
    for i=1:length(tobs)
      xmod[i],ymod[i] = return_xy(tobs[i],elements)
    end
    chisq[j,k] =sum((xobs .- xmod).^2+(yobs .- ymod).^2)
    if chisq[j,k] < chi_best
      chi_best = chisq[j,k]
      ibest = j
      elements_best = elements
      println(j," chi_best: ",chi_best," elements_best: ",elements_best)
      clf()
      plot(xobs,xmod)
      plot(yobs,ymod)
    end
  end
end


xmod, ymod , res = fit_keplerian(tobs,xobs,yobs,elements_best)
clf()
#plot(tobs,xobs-xmod)
#plot(tobs,yobs-ymod)

# Plot the difference in angle between the n-body and Keplerian angles:
dtheta = atan2(yobs.*xmod .- ymod.*xobs,xobs.*xmod+yobs.*ymod)
plot(tobs,dtheta)
