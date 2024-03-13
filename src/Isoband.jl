module Isoband

using isoband_jll

export isobands, isolines

struct ReturnValue{T}
    xs::Ptr{T}
    ys::Ptr{T}
    ids::Ptr{Cint}
    len::Cint
end

isoband_float_type(_::Type{Float64}) = Cdouble
isoband_float_type(_::Type{Float32}) = Cfloat

float_eltype(x) = eltype(x) |> float

"""
    isobands(xs, ys, zs, low::Real, high::Real)

Create one isoband from the matrix `zs` for the boundaries `low` and `high`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Returns a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one polygon, the polygons can be outer polygons or
holes and are given in no particular order. Therefore, they must probably be post-processed to
feed them to plotting packages.
"""
isobands(xs, ys, zs, low::Real, high::Real) = isobands(xs, ys, zs, float(low), float(high))

function isobands(xs, ys, zs, low::T, high::T) where {T <: AbstractFloat}
    results = isobands(xs, ys, zs, T[low], T[high])
    results[1]
end

function isobands(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix, lows::AbstractVector, highs::AbstractVector)
    F = promote_type(float_eltype(xs), float_eltype(ys), float_eltype(zs), float_eltype(lows), float_eltype(highs))
    isobands(F.(xs), F.(ys), F.(zs), F.(lows), F.(highs))
end

"""
    isobands(xs::Vector{Float64}, ys::Vector{Float64}, zs::Matrix{Float64}, low_values::Vector{Float64}, high_values::Vector{Float64})

Create a vector of isobands from the matrix `zs` for all pairs in `low_values` and `high_values`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Each entry of the return vector is a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one polygon, the polygons can be outer polygons or
holes and are given in no particular order. Therefore, they must probably be post-processed to
feed them to plotting packages.
"""
function isobands(
    xs::Vector{T},
    ys::Vector{T},
    zs::Matrix{T},
    low_values::Vector{T},
    high_values::Vector{T},
) where {T <: AbstractFloat}
    lenx = length(xs)
    leny = length(ys)
    nrow, ncol = size(zs)

    lenx != ncol && error("Length of x ($(length(xs))) must be equal to number of columns in z ($(size(zs, 2)))")
    leny != nrow && error("Length of y $(length(ys)) must be equal to number of rows in z ($(size(zs, 1))")

    length(low_values) != length(high_values) && error("Number of low values ($(length(low_values)))and high values ($(length(high_values))) must be equal.")

    nbands = length(low_values)

    result = if (F = isoband_float_type(T)) ≡ Float64
        ccall((:isobands_impl, libisoband), Ptr{ReturnValue{Cdouble}}, (
                Ptr{Cdouble},
                Cint,
                Ptr{Cdouble},
                Cint,
                Ptr{Cdouble},
                Cint,
                Cint,
                Ptr{Cdouble},
                Ptr{Cdouble},
                Cint,
            ),
            xs, length(xs),
            ys, length(ys),
            zs, size(zs, 1), size(zs, 2),
            low_values, high_values, nbands
        )
    elseif F ≡ Float32
        ccall((:isobands32_impl, libisoband), Ptr{ReturnValue{Cfloat}}, (
                Ptr{Cfloat},
                Cint,
                Ptr{Cfloat},
                Cint,
                Ptr{Cfloat},
                Cint,
                Cint,
                Ptr{Cfloat},
                Ptr{Cfloat},
                Cint,
            ),
            xs, length(xs),
            ys, length(ys),
            zs, size(zs, 1), size(zs, 2),
            low_values, high_values, nbands
        )
    else
        throw(ArgumentError("$F is invalid"))
    end
    
    returnvalues = unsafe_wrap(Vector{ReturnValue{F}}, result, nbands, own = true)

    groups = map(returnvalues) do rv
        n = rv.len
        xsr = unsafe_wrap(Vector{F}, rv.xs, n, own = true)
        ysr = unsafe_wrap(Vector{F}, rv.ys, n, own = true)
        idr = unsafe_wrap(Vector{Cint}, rv.ids, n, own = true)
        (x = xsr, y = ysr, id = idr)
    end

    groups
end

"""
    isolines(xs, ys, zs, value::Real)

Create one isoline from the matrix `zs` for `value`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Returns a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one line.
"""

isolines(xs, ys, zs, value::Real) = isolines(xs, ys, zs, float(value))

function isolines(xs, ys, zs, value::T) where {T <: AbstractFloat}
    results = isolines(xs, ys, zs, T[value])
    results[1]
end

function isolines(xs::AbstractVector, ys::AbstractVector, zs::AbstractMatrix, values::AbstractVector)
    F = promote_type(float_eltype(xs), float_eltype(ys), float_eltype(zs), float_eltype(values))
    isolines(F.(xs), F.(ys), F.(zs), F.(values))
end

"""
    isolines(xs::Vector{Float64}, ys::Vector{Float64}, zs::Matrix{Float64}, values::Vector{Float64})

Create a vector of isolines from the matrix `zs` for all `values`.
The rows of `zs` correspond to the linear spaced values in `ys` and the columns to `xs`.
Each entry of the return vector is a NamedTuple with two Vector{Float64} fields `x` and `y`, and the
Vector{Int} field `id`. Each unique id marks one line.
"""
function isolines(
    xs::Vector{T},
    ys::Vector{T},
    zs::Matrix{T},
    values::Vector{T},
) where {T <: AbstractFloat}
    lenx = length(xs)
    leny = length(ys)
    nrow, ncol = size(zs)

    lenx != ncol && error("Length of x ($(length(xs))) must be equal to number of columns in z ($(size(zs, 2)))")
    leny != nrow && error("Length of y $(length(ys)) must be equal to number of rows in z ($(size(zs, 1))")

    nvalues = length(values)

    result = if (F = isoband_float_type(T)) ≡ Float64
        ccall((:isolines_impl, libisoband), Ptr{ReturnValue{Cdouble}}, (
                Ptr{Cdouble},
                Cint,
                Ptr{Cdouble},
                Cint,
                Ptr{Cdouble},
                Cint,
                Cint,
                Ptr{Cdouble},
                Cint,
            ),
            xs, length(xs),
            ys, length(ys),
            zs, size(zs, 1), size(zs, 2),
            values, nvalues
        )
    elseif F ≡ Float32
        ccall((:isolines32_impl, libisoband), Ptr{ReturnValue{Cfloat}}, (
                Ptr{Cfloat},
                Cint,
                Ptr{Cfloat},
                Cint,
                Ptr{Cfloat},
                Cint,
                Cint,
                Ptr{Cfloat},
                Cint,
            ),
            xs, length(xs),
            ys, length(ys),
            zs, size(zs, 1), size(zs, 2),
            values, nvalues
        )
    else
        throw(ArgumentError("$F is invalid"))
    end

    returnvalues = unsafe_wrap(Vector{ReturnValue{F}}, result, nvalues, own = true)

    groups = map(returnvalues) do rv
        n = rv.len
        xsr = unsafe_wrap(Vector{F}, rv.xs, n, own = true)
        ysr = unsafe_wrap(Vector{F}, rv.ys, n, own = true)
        idr = unsafe_wrap(Vector{Cint}, rv.ids, n, own = true)
        (x = xsr, y = ysr, id = idr)
    end
    
    groups
end


end
