# This code tests two functions: keplerij! and kep_elliptic!

@testset "keplerij" begin

# Next, try computing two-body Keplerian Jacobian:

t0 = 7257.93115525
elements = "elements.txt"
system = [3,1,1]

n = 3
h  = 0.05
hbig  = big(h)
tmax = 600.0
#dlnq = 1e-8
dlnq = 1e-20
init = IC(elements,system)
#elements[2,1] = 0.75
init.elements[2,1] = 1.0
init.elements[3,1] = 1.0

x0=zeros(NDIM,n)
v0=zeros(NDIM,n)

for iter = 1:2
@time x0,v0 = init_nbody(init,t0)
h = 0.05; hbig = big(h)
 if iter == 2
   # Reduce masses to trigger hyperbolic routine:
    init.m[1:n] *= 1e-1
    h = 0.05
    hbig = big(h)
 end
# Tilt the orbits a bit:
x0[2,1] = 5e-1*sqrt(x0[1,1]^2+x0[3,1]^2)
x0[2,2] = -5e-1*sqrt(x0[1,2]^2+x0[3,2]^2)
v0[2,1] = 5e-1*sqrt(v0[1,1]^2+v0[3,1]^2)
v0[2,2] = -5e-1*sqrt(v0[1,2]^2+v0[3,2]^2)

jac_ij = zeros(Float64,14,14)
dqdt_ij = zeros(Float64,14)
dqdt_num = zeros(BigFloat,14)
i=1 ; j=2
x = copy(x0) ; v=copy(v0)
# Predict values of s:
#println("Before first step: ",x," ",v)
keplerij!(init.m,x,v,i,j,h,jac_ij,dqdt_ij)
#println("After first step: ",x," ",v)
x0 = copy(x) ; v0 = copy(v)
xbig = big.(x) ; vbig=big.(v); mbig = big.(init.m)
keplerij!(init.m,x,v,i,j,h,jac_ij,dqdt_ij)
# Now compute Jacobian with BigFloat precision:
jac_ij_big = zeros(BigFloat,14,14)
dqdt_ij_big = zeros(BigFloat,14)
KEPLER_TOL = sqrt(eps(big(1.0)))
keplerij!(mbig,xbig,vbig,i,j,hbig,jac_ij_big,dqdt_ij_big)
#println("jac_ij: ",convert(Array{Float64,2},jac_ij_big))
#println("jac_ij - jac_ij_big: ",convert(Array{Float64,2},jac_ij_big)-jac_ij)
println("max(jac_ij - jac_ij_big): ",maxabs(convert(Array{Float64,2},jac_ij_big)-jac_ij))

# xtest = copy(x0) ; vtest=copy(v0)
# keplerij!(m,xtest,vtest,i,j,h,jac_ij)
# println("Test of jacobian vs. none: ",maximum(abs(x-xtest)),maximum(abs(v-vtest)))

# Now, compute the derivatives numerically:
jac_ij_num = zeros(BigFloat,14,14)
global xsave = big.(x)
vsave = big.(v)
msave = big.(init.m)

# Compute the time derivatives:
# Initial positions, velocities & masses:
xm = big.(x0)
vm = big.(v0)
mm = big.(msave)
dq = dlnq * hbig
hbig -= dq
keplerij!(mm,xm,vm,i,j,hbig)
xp = big.(x0)
vp = big.(v0)
hbig += 2dq
keplerij!(mm,xp,vp,i,j,hbig)
# Now x & v are final positions & velocities after time step
for k=1:3
  dqdt_num[   k] = .5*(xp[k,i]-xm[k,i])/dq
  dqdt_num[ 3+k] = .5*(vp[k,i]-vm[k,i])/dq
  dqdt_num[ 7+k] = .5*(xp[k,j]-xm[k,j])/dq
  dqdt_num[10+k] = .5*(vp[k,j]-vm[k,j])/dq
end
hbig = big(h)

for jj=1:3
  # Initial positions, velocities & masses:
  xm = big.(x0)
  vm = big.(v0)
  mm = big.(msave)
  dq = dlnq * xm[jj,i]
  if xm[jj,i] != 0.0
    xm[jj,i] -=  dq
  else
    dq = dlnq
    xm[jj,i] = -dq
  end
  keplerij!(mm,xm,vm,i,j,hbig)
  xp = big.(x0)
  vp = big.(v0)
  if xm[jj,i] != 0.0
    xp[jj,i] +=  dq
  else
    dq = dlnq
    xp[jj,i] = dq
  end
  keplerij!(mm,xp,vp,i,j,hbig)
  # Now x & v are final positions & velocities after time step
  for k=1:3
    jac_ij_num[   k,  jj] = .5*(xp[k,i]-xm[k,i])/dq
    jac_ij_num[ 3+k,  jj] = .5*(vp[k,i]-vm[k,i])/dq
    jac_ij_num[ 7+k,  jj] = .5*(xp[k,j]-xm[k,j])/dq
    jac_ij_num[10+k,  jj] = .5*(vp[k,j]-vm[k,j])/dq
  end
  xm = big.(x0)
  vm = big.(v0)
  mm  = big.(msave)
  dq = dlnq * vm[jj,i]
  if vm[jj,i] != 0.0
    vm[jj,i] -=  dq
  else
    dq = dlnq
    vm[jj,i] = -dq
  end
  keplerij!(mm,xm,vm,i,j,hbig)
  xp = big.(x0)
  vp = big.(v0)
  mm  = big.(msave)
  if vp[jj,i] != 0.0
    vp[jj,i] +=  dq
  else
    dq = dlnq
    vp[jj,i] = dq
  end
  keplerij!(mm,xp,vp,i,j,hbig)
  for k=1:3
    jac_ij_num[   k,3+jj] = .5*(xp[k,i]-xm[k,i])/dq
    jac_ij_num[ 3+k,3+jj] = .5*(vp[k,i]-vm[k,i])/dq
    jac_ij_num[ 7+k,3+jj] = .5*(xp[k,j]-xm[k,j])/dq
    jac_ij_num[10+k,3+jj] = .5*(vp[k,j]-vm[k,j])/dq
  end
end

# Now vary mass of inner planet:
xm= big.(x0)
vm= big.(v0)
mm= big.(msave)
dq = mm[i]*dlnq
mm[i] -= dq
keplerij!(mm,xm,vm,i,j,hbig)
xp= big.(x0)
vp= big.(v0)
mp= big.(msave)
dq = mp[i]*dlnq
mp[i] += dq
keplerij!(mp,xp,vp,i,j,hbig)
for k=1:3
  jac_ij_num[   k,7] = .5*(xp[k,i]-xm[k,i])/dq
  jac_ij_num[ 3+k,7] = .5*(vp[k,i]-vm[k,i])/dq
  jac_ij_num[ 7+k,7] = .5*(xp[k,j]-xm[k,j])/dq
  jac_ij_num[10+k,7] = .5*(vp[k,j]-vm[k,j])/dq
end
# The mass doesn't change:
jac_ij_num[7,7] =  1.0

for jj=1:3
  # Now vary parameters of outer planet:
  xm = big.(x0)
  vm = big.(v0)
  mm = big.(msave)
  dq = dlnq * xm[jj,j]
  if xm[jj,j] != 0.0
    xm[jj,j] -=  dq
  else
    dq = dlnq
    xm[jj,j] = -dq
  end
  keplerij!(mm,xm,vm,i,j,hbig)
  xp = big.(x0)
  vp = big.(v0)
  if xp[jj,j] != 0.0
    xp[jj,j] +=  dq
  else
    dq = dlnq
    xp[jj,j] = dq
  end
  keplerij!(mm,xp,vp,i,j,hbig)
  for k=1:3
    jac_ij_num[   k,7+jj] = .5*(xp[k,i]-xm[k,i])/dq
    jac_ij_num[ 3+k,7+jj] = .5*(vp[k,i]-vm[k,i])/dq
    jac_ij_num[ 7+k,7+jj] = .5*(xp[k,j]-xm[k,j])/dq
    jac_ij_num[10+k,7+jj] = .5*(vp[k,j]-vm[k,j])/dq
  end
  xm= big.(x0)
  vm= big.(v0)
  mm = big.(msave)
  dq = dlnq * vm[jj,j]
  if vm[jj,j] != 0.0
    vm[jj,j] -=  dq
  else
    dq = dlnq
    vm[jj,j] = -dq
  end
  keplerij!(mm,xm,vm,i,j,hbig)
  xp= big.(x0)
  vp= big.(v0)
  if vp[jj,j] != 0.0
    vp[jj,j] +=  dq
  else
    dq = dlnq
    vp[jj,j] = dq
  end
  keplerij!(mm,xp,vp,i,j,hbig)
  for k=1:3
    jac_ij_num[   k,10+jj] = .5*(xp[k,i]-xm[k,i])/dq
    jac_ij_num[ 3+k,10+jj] = .5*(vp[k,i]-vm[k,i])/dq
    jac_ij_num[ 7+k,10+jj] = .5*(xp[k,j]-xm[k,j])/dq
    jac_ij_num[10+k,10+jj] = .5*(vp[k,j]-vm[k,j])/dq
  end
end

# Now vary mass of outer planet:
xm = big.(x0)
vm = big.(v0)
mm = big.(msave)
dq = mm[j]*dlnq
mm[j] -= dq
keplerij!(mm,xm,vm,i,j,hbig)
xp = big.(x0)
vp = big.(v0)
mp = big.(msave)
dq = mp[j]*dlnq
mp[j] += dq
keplerij!(mp,xp,vp,i,j,hbig)
for k=1:3
  jac_ij_num[   k,14] = .5*(xp[k,i]-xm[k,i])/dq
  jac_ij_num[ 3+k,14] = .5*(vp[k,i]-vm[k,i])/dq
  jac_ij_num[ 7+k,14] = .5*(xp[k,j]-xm[k,j])/dq
  jac_ij_num[10+k,14] = .5*(vp[k,j]-vm[k,j])/dq
end
# The mass doesn't change:
jac_ij_num[14,14] =  1.0

#println(jac_ij)
#println(jac_ij_num)
#println(jac_ij./jac_ij_num)
global emax = 0.0; global imax = 0; global jmax = 0
for i=1:14, j=1:14
  if jac_ij[i,j] != 0.0
    diff = abs(convert(Float64,jac_ij_num[i,j])/jac_ij[i,j]-1.0)
    if  diff > emax
      emax = diff; imax = i; jmax = j
    end
  end
end
println("Maximum fractional error: ",emax," ",imax," ",jmax)
#println(jac_ij)
#println(convert(Array{Float64,2},jac_ij_num))
println("Maximum jac_ij error:   ",maxabs(convert(Array{Float64,2},jac_ij_num)-jac_ij))
println("Maximum jac_ij_big-jac_ij_num:   ",maxabs(convert(Array{Float64,2},jac_ij_num-jac_ij_big)))
println("Max dqdt error: ",maxabs(dqdt_ij-convert(Array{Float64,1},dqdt_num)))

#@test isapprox(jac_ij_num,jac_ij;norm=maxabs)
@test isapprox(dqdt_ij,convert(Array{Float64,1},dqdt_num);norm=maxabs)
#println("dqdt: ",dqdt_ij," ",dqdt_num," Max error: ",maximum(dqdt_ij-dqdt_num))
end
end
