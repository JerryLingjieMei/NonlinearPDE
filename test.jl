include("./src/solver.jl")
using PyPlot

#=
solver.jl:
- Julia version: 0.6.2
- Author: jerrymei
- Date: 2018-04-27
=#

k(x, y)=0.6 * sin(2 * pi * x) + 0.2 * sin(4 * pi * y) + 0.1 * sin(10 * pi * x + 8 * pi * y)
N_mesh=512

x=linspace(0,1 - 1 / N_mesh,N_mesh)
y=linspace(0,1 - 1 / N_mesh,N_mesh)
sol=zeros(N_mesh,N_mesh)
for i=1:N_mesh
    for j=1:N_mesh
        sol[i,j]=k(x[i],y[j])
    end
end

f=get_f(sol,1.)

u,err=solver(f,1.,1e-6)
figure(1)
plot_surface(linspace(0,1 - 1 / N_mesh,N_mesh),linspace(0,1 - 1 / N_mesh,N_mesh),u)
println("Norm of the distance between Du and f ", vecnorm(get_f(u,1.) - f)/vecnorm(f))

u,err=multigrid_solver(f,1.,1e-6)
figure(1)
plot_surface(linspace(0,1 - 1 / N_mesh,N_mesh),linspace(0,1 - 1 / N_mesh,N_mesh),u)
println("Norm of the distance between Du and f ", vecnorm(get_f(u,1.) - f)/vecnorm(f))