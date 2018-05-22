#=
multigrid:
- Julia version: 0.6.2
- Author: jerrymei
- Date: 2018-05-22
=#

mutable struct gridWrapper{T}
    u::Array{T,2}
    newu::Array{T,2}
    dx::Array{T,2}
    dy::Array{T,2}
    rhs::Array{T,2}
    f::Array{T,2}
    fft::rFFTWrapper
end

function gridWrapper{T}(f::AbstractArray{T,2})
    fft=rFFTWrapper(f)
    u=zeros(f)
    newu=similar(u)
    dx=similar(u)
    dy=similar(u)
    rhs=similar(u)
    return gridWrapper(u,newu,dx,dy,rhs,f,fft)
end


function restriction!{T}(fine::AbstractArray{T,2}, coarse::AbstractArray{T,2})
    nx,ny=size(coarse)
    @inbounds for x in 1:nx
        mx=2 * x - 1
        lx=mx > 1?mx - 1:2 * nx
        ux=mx + 1
        for y in 1:ny
            my=2 * y - 1
            ly=my > 1?my - 1:2 * ny
            uy=my + 1
            coarse[x,y]=(fine[lx,ly] + 2. * fine[mx,ly] + fine[ux,ly])/16.+ (fine[lx,my] + 2. * fine[mx,my] + fine[ux,my])/8.+ (fine[lx,uy] + 2. * fine[mx,uy] + fine[ux,uy]) / 16.
        end
    end
end

function interpolation!{T}(coarse::AbstractArray{T,2}, fine::AbstractArray{T,2})
    nx,ny =size(coarse)
    @inbounds for x in 1:nx
        mx=2 * x - 1
        ux=x < nx?x + 1:1
        for y in 1:ny
            my=2 * y - 1
            uy=y < ny?y + 1:1
            fine[mx,my]=coarse[x,y]
            fine[mx,my + 1]=(coarse[x,y] + coarse[x,uy]) / 2.
            fine[mx + 1,my]=(coarse[x,y] + coarse[ux,y]) / 2.
            fine[mx + 1,my + 1]=(coarse[x,y] + coarse[x,uy] + coarse[ux,y] + coarse[ux,uy]) / 4.
        end
    end
end


function full_cycle{T}(iteration,grid1,grid2,grid3,grid4,eps::T,tol::T)
    iteration(grid4,eps)
    interpolation!(grid4.u,grid3.u)
    iteration(grid3,eps)
    restriction!(grid3.u,grid4.u)

    iteration(grid4,eps)
    interpolation!(grid4.u,grid3.u)
    iteration(grid3,eps)
    interpolation!(grid3.u,grid2.u)
    iteration(grid2,eps)
    restriction!(grid2.u,grid3.u)
    iteration(grid3,eps)
    restriction!(grid3.u,grid4.u)

    iteration(grid4,eps)
    interpolation!(grid4.u,grid3.u)
    iteration(grid3,eps)
    interpolation!(grid3.u,grid2.u)
    iteration(grid2,eps)
    interpolation!(grid2.u,grid1.u)
end