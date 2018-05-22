#=
fft:
- Julia version: 0.6.2
- Author: jerrymei
- Date: 2018-05-05
=#
using FFTW

mutable struct rFFTWrapper{T}
    src::Array{T,2}
    src_hat::Array{Complex{T},2}
    dest_hat::Array{Complex{T},2}
    dest::Array{T,2}
    plan::FFTW.Plan
    inv_plan::FFTW.Plan
    src_hat_calculated::Bool
end

function rFFTWrapper{T}(input::AbstractArray{T,2})
    mx,my=size(input)
    nx,ny=div(mx,2) + 1,my
    src=input
    src_hat=Array{Complex{T}}(nx,ny)
    dest_hat=similar(src_hat)
    dest=similar(input)
    plan=plan_rfft(src)
    inv_plan=plan_irfft(src_hat,mx)
    return rFFTWrapper(src,src_hat,dest_hat,dest,plan,inv_plan,false)
end

function set_src!{T}(fft::rFFTWrapper, input::AbstractArray{T,2})
    if !(input === fft.src)
        assert(size(input) == size(fft.src))
        fft.src=input
        fft.src_hat_calculated=false
    end
end

function set_dest!{T}(fft::rFFTWrapper, output::AbstractArray{T,2})
    assert(size(output) == size(fft.dest))
    fft.dest=output
end

function set_not_calculated(fft::rFFTWrapper)
    fft.src_hat_calculated=false
end

function derivx!(fft::rFFTWrapper)
    mx,my=size(fft.src)
    nx,ny=size(fft.src_hat)
    if !fft.src_hat_calculated
        A_mul_B!(fft.src_hat,fft.plan,fft.src)
        fft.src_hat_calculated=true
    end
    @inbounds for x in 0:nx - 1
        for y in 0:ny - 1
            fft.dest_hat[x + 1,y + 1]=fft.src_hat[x + 1,y + 1] * x * 2pi * im
        end
    end
    A_mul_B!(fft.dest,fft.inv_plan,fft.dest_hat)
end

function derivy!(fft::rFFTWrapper)
    mx,my=size(fft.src)
    nx,ny=size(fft.src_hat)
    if !fft.src_hat_calculated
        A_mul_B!(fft.src_hat,fft.plan,fft.src)
        fft.src_hat_calculated=true
    end
    @inbounds for x in 0:nx - 1
        for y in 0:ny - 1
            fft.dest_hat[x + 1,y + 1]= fft.src_hat[x + 1,y + 1] * (2 * y < my?y:y - my) * 2pi * im
        end
    end
    A_mul_B!(fft.dest,fft.inv_plan,fft.dest_hat)
end

function laplacian!(fft::rFFTWrapper)
    mx,my=size(fft.src)
    nx,ny=size(fft.src_hat)
    if !fft.src_hat_calculated
        A_mul_B!(fft.src_hat,fft.plan,fft.src)
        fft.src_hat_calculated=true
    end
    @inbounds for x in 0:nx - 1
        for y in 0:ny - 1
            fft.dest_hat[x + 1,y + 1]= fft.src_hat[x + 1,y + 1] * (x^2 + (2 * y < my?y:y - my)^2) * (2pi * im)^2
        end
    end
    A_mul_B!(fft.dest,fft.inv_plan,fft.dest_hat)
end

function solve_laplacian!(fft::rFFTWrapper)
    mx,my=size(fft.src)
    nx,ny=size(fft.src_hat)
    if !fft.src_hat_calculated
        A_mul_B!(fft.src_hat,fft.plan,fft.src)
        fft.src_hat_calculated=true
    end
    @inbounds for x in 0:nx - 1
        for y in 1:ny - 1
            fft.dest_hat[x + 1,y + 1]= fft.src_hat[x + 1,y + 1] / (x^2 + (2 * y < my?y:y - my)^2) / (2pi * im)^2
        end
    end
    @inbounds for x in 1:nx - 1
        fft.dest_hat[x + 1,1]=fft.src_hat[x + 1,1] / x^2 / (2pi * im)^2
    end
    A_mul_B!(fft.dest,fft.inv_plan,fft.dest_hat)
end

function derivx!{T}(input::AbstractArray{T,2}, output::AbstractArray{T,2}, fft::rFFTWrapper)
"""Compute y derivative"""
    set_src!(fft,input)
    set_dest!(fft,output)
    derivx!(fft)
end

function derivy!{T}(input::AbstractArray{T,2}, output::AbstractArray{T,2}, fft::rFFTWrapper)
"""Compute y derivative"""
    set_src!(fft,input)
    set_dest!(fft,output)
    derivy!(fft)
end

function laplacian!{T}(input::AbstractArray{T,2}, output::AbstractArray{T,2}, fft::rFFTWrapper)
"""Compute the Laplacian"""
    set_src!(fft,input)
    set_dest!(fft,output)
    laplacian!(fft)
end

function solve_laplacian!{T}(input::AbstractArray{T,2}, output::AbstractArray{T,2}, fft::rFFTWrapper)
"""Solve Poisson equation"""
    set_src!(fft,input)
    set_dest!(fft,output)
    solve_laplacian!(fft)
end


