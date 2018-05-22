include("fft.jl")
include("multigrid.jl")
function get_f{T}(u::AbstractArray{T,2}, eps::T)
"""Get the f with left hand side"""
    mx,my=size(u)
    fft=rFFTWrapper(u)
    lap = similar(u)
    laplacian!(u,lap,fft)
    dx=similar(u)
    dy=similar(u)
    set_not_calculated(fft)
    derivx!(u,dx,fft)
    derivy!(u,dy,fft)
    @inbounds for x in 1:mx
        for y in 1:my
            lap[x,y] += eps * (dx[x,y]^2 + dy[x,y]^2 + lap[x,y] * u[x,y])
        end
    end
    return lap
end

function solver{T}(f::AbstractArray{T,2}, eps::T, tol::T, maxiter=50)
"""Solving the equation with laplacian"""
    mx,my=size(f)
    fft=rFFTWrapper(f)
    u=zeros(f)
    newu=similar(u)
    dx=similar(u)
    dy=similar(u)
    rhs=similar(u)
    iter=0
    err=T[]
    while iter < maxiter
        iter += 1
        derivx!(u,dx,fft)
        derivy!(u,dy,fft)
        @inbounds for x in 1:mx
            for y in 1:my
                rhs[x,y]= (f[x,y] - eps * (dx[x,y]^2 + dy[x,y]^2)) / (one(T) + eps * u[x,y])
            end
        end
        solve_laplacian!(rhs,newu,fft)
        e=vecnorm(newu - u) / vecnorm(f)
        push!(err,e)
        if e < tol
            break
        end
        u,newu=newu,u
    end
    return u,err
end

function iteration{T}(grid::gridWrapper, eps::T, err::Array{T,1}, maxiter=4)
    mx,my=size(grid.f)
    iter=0
    while iter < maxiter
        iter += 1
        derivx!(grid.u,grid.dx,grid.fft)
        derivy!(grid.u,grid.dy,grid.fft)
        @inbounds for x in 1:mx
            for y in 1:my
                grid.rhs[x,y]= (grid.f[x,y] - eps * (grid.dx[x,y]^2 + grid.dy[x,y]^2)) / (one(T) + eps * grid.u[x,y])
            end
        end
        solve_laplacian!(grid.rhs,grid.newu,grid.fft)
        grid.u,grid.newu=grid.newu,grid.u
        push!(err,vecnorm(grid.newu - grid.u) / vecnorm(grid.f))
    end
end

function iteration{T}(grid::gridWrapper, eps::T, maxiter=4)
    mx,my=size(grid.f)
    iter=0
    while iter < maxiter
        iter += 1
        derivx!(grid.u,grid.dx,grid.fft)
        derivy!(grid.u,grid.dy,grid.fft)
        @inbounds for x in 1:mx
            for y in 1:my
                grid.rhs[x,y]= (grid.f[x,y] - eps * (grid.dx[x,y]^2 + grid.dy[x,y]^2)) / (one(T) + eps * grid.u[x,y])
            end
        end
        solve_laplacian!(grid.rhs,grid.newu,grid.fft)
        grid.u,grid.newu=grid.newu,grid.u
    end
end


function multigrid_solver{T}(f::AbstractArray{T,2}, eps::T, tol::T, maxiter=32)
"""Solve the equation with multigrid"""
    mx,my=size(f)
    grid1=gridWrapper(f)
    f2=similar(f,div(mx,2),div(my,2))
    restriction!(f,f2)
    grid2=gridWrapper(f2)
    f3=similar(f,div(mx,4),div(my,4))
    restriction!(f2,f3)
    grid3=gridWrapper(f3)
    f4=similar(f,div(mx,8),div(my,8))
    restriction!(f3,f4)
    grid4=gridWrapper(f4)
    err=T[]
    full_cycle(iteration,grid1,grid2,grid3,grid4,eps,tol)
    iter=0
    while iter < maxiter
        iter += 1
        iteration(grid1,eps,err)
        if vecnorm(grid1.u - grid1.newu) < tol * vecnorm(grid1.f)
            break
        end
    end
    return grid1.u,err
end

