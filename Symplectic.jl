####################################
##  Types for symplectic vectors  ##
####################################

abstract type AbstractSymplecticVector end

struct SymplecticVector{n, d} <: AbstractSymplecticVector
    z::FqMatrix
    x::FqMatrix
    function SymplecticVector{n, d}(z::FqMatrix, x::FqMatrix) where {n, d}
        new(z, x)
    end
end

# Outer constructor accepting integers and converting them to finite field elements
function SymplecticVector{n, d}(z::Array, x::Array) where {n, d}
    FF = finite_field(d)[1]
    MM = matrix_space(FF, n, 1)
    return SymplecticVector{n, d}(MM(z), MM(x))
end

# Overloading builtin functions for SymplecticVector types
function hash(v::SymplecticVector{n, d}) where {n, d}
    return hash((v.z, v.x))
end

# Printing SymplecticVectors
function show(io::IO, v::SymplecticVector{n, d}) where {n, d}
    print(io, reduce(*, [
            "[",
            join(permutedims(hcat(string.(v.z[:,1]), repeat([" "], n)))[:])[1:end-1],
            "|",
            join(permutedims(hcat(string.(v.x[:,1]), repeat([" "], n)))[:])[1:end-1],
            "]"
        ]))
end

# Functions for generating random symplectic vectors
function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, d}}) where {n, d}
    FF = finite_field(d)[1]
    return SymplecticVector{n, d}(rand(FF, n), rand(FF, n))
end

function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, d}}, dims...) where {n, d}
    FF = finite_field(d)[1]
    return SymplecticVector{n, d}.(rand(FF, n, dims...), rand(FF, n, dims...))
end

# Functions for extracting information from SymplecticVector types
dimension(v::SymplecticVector{n, d}) where {n, d} = 2n
halfdimension(v::SymplecticVector{n, d}) where {n, d} = n

# Basic arithmetic with symplectic vectors
==(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = (u.z == v.z) && (u.x == v.x)

+(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = SymplecticVector{n, d}(u.z+v.z, u.x+v.x)
-(v::SymplecticVector{n, d}) where {n, d} = SymplecticVector{n, d}(-v.z, -v.x)
-(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = u + (-v)
*(k::FqFieldElem, v::SymplecticVector{n, d}) where {n, d} = SymplecticVector{n, d}(k * v.z, k * v.x)
*(v::SymplecticVector{n, d}, k::FqFieldElem) where {n, d} = SymplecticVector{n, d}(k * v.z, k * v.x)

function innerproduct(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    (transpose(u.z) * v.z + transpose(u.x) * v.x)[1]
end
⋅ = innerproduct

function symplecticform(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    (transpose(u.z) * v.x - transpose(u.x) * v.z)[1]
end
⋆ = symplecticform