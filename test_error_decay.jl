using JLD
include("src/solver.jl")
k(x, y)=0.6 * sin(2 * pi * x) + 0.2 * sin(4 * pi * y) + 0.1 * sin(10 * pi * x + 8 * pi * y)
N_mesh=512
arr1=[]
arr2=[]
for N_mesh in [256, 512, 1024, 2048, 4096]
    println(N_mesh)
    x=linspace(0,1 - 1 / N_mesh,N_mesh)
    y=linspace(0,1 - 1 / N_mesh,N_mesh)
    sol=zeros(N_mesh,N_mesh)  
    for i=1:N_mesh
        for j=1:N_mesh
            sol[i,j]=k(x[i],y[j])
        end
    end
    f=get_f(sol,1.)
    u,err1=solver(f,1.,1e-6)
    push!(arr1,err1)
    u,err2=multigrid_solver(f,1.,1e-6)
    push!(arr2,err2)
end
brr1=[]
brr2=[]

for N_mesh in [256, 512, 1024, 2048, 4096]
    println(N_mesh)
    x=linspace(0,1 - 1 / N_mesh,N_mesh)
    y=linspace(0,1 - 1 / N_mesh,N_mesh)
    sol=zeros(N_mesh,N_mesh)  
    for i=1:N_mesh
        for j=1:N_mesh
            sol[i,j]=k(x[i],y[j])
        end
    end
    f=get_f(sol,1.)
    push!(brr1,@elapsed solver(f,1.,1e-6))
    push!(brr2,@elapsed multigrid_solver(f,1.,1e-6))
end

save("plotdata.jld","arr1",arr1,"arr2",arr2,"brr1",brr1,"brr2",brr2)