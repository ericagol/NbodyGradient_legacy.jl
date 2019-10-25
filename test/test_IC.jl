# Testing outputs of functions with different precisions
@testset "Test IC" begin

elements = "elements.txt"
system = [4,1,1,1]

t0 = 7257.93115525
t0big = big(t0)
h  = 0.04
hbig = big(0.04)
tmax = 100.0
tmaxbig = big(1000.0)
rstar = 1e12
rstarbig = big(1e12)

init = IC(elements,system)
initbig = IC(elements,system;prec=BigFloat)

@test isapprox(init.elements,initbig.elements)

x,v,jac_init = init_nbody(init,t0)
xbig,vbig,jac_init_big = init_nbody(initbig,t0big)

@test isapprox(x,xbig)
@test isapprox(v,vbig)
@test isapprox(jac_init,jac_init_big)

ntt = zeros(Int64,init.nbody)
for i=2:init.nbody
  ntt[i] = ceil(Int64,tmax/init.elements[i,2])+3
end
tt  = zeros(init.nbody,maximum(ntt))
count = zeros(Int64,init.nbody)
dq = NbodyGradient.ttv_elements!(init,t0,h,tmax,tt,count,0.0,0,0,rstar)

ntt = zeros(Int64,init.nbody)
for i=2:init.nbody
  ntt[i] = ceil(Int64,tmax/initbig.elements[i,2])+3
end
ttbig  = zeros(BigFloat,init.nbody,maximum(ntt))
countbig = zeros(Int64,init.nbody)
dqbig = NbodyGradient.ttv_elements!(initbig,t0big,hbig,tmaxbig,ttbig,countbig,big(0.0),0,0,rstarbig)

@test isapprox(dqbig,dq)

end
