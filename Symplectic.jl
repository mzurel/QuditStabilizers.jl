################################################
##  Types and methods for symplectic vectors  ##
################################################

abstract type AbstractSymplecticVector end

"""
    SymplecticVector{n, d}(z::FqMatrix, x::FqMatrix)

Construct a symplectic vector over ``ℤₚ`` (p=d) with `z` and `x` components as `FqMatrix` objects.

# Arguments
- `z::FqMatrix`: The first component of the symplectic vector.
- `x::FqMatrix`: The second component of the symplectic vector.

# Example
```jldoctest
julia> FF = finite_field(3)[1]
Prime field of characteristic 3

julia> MM = matrix_space(FF, 2, 1)
Matrix space of 2 rows and 1 column
  over prime field of characteristic 3

julia> z = MM([1;2])
[1]
[2]

julia> x = MM([2;0])
[2]
[0]

julia> v = SymplecticVector{2, 3}(z, x)
[1 2|2 0]
```
"""
struct SymplecticVector{n, d} <: AbstractSymplecticVector
    z::FqMatrix
    x::FqMatrix
    function SymplecticVector{n, d}(z::FqMatrix, x::FqMatrix) where {n, d}
        new(z, x)
    end
end

"""
    SymplecticVector{n, d}(z::Array, x::Array)

Outer constructor for a symplectic vector with `z` and `x` components provided as `Array`
objects, automatically converting them to `FqMatrix` objects over the finite field `F_q`.

# Example
```jldoctest
julia> v = SymplecticVector{2, 3}([1,2], [2,0])
[1 2|2 0]
```
"""
function SymplecticVector{n, d}(z::Array, x::Array) where {n, d}
    FF = finite_field(d)[1]
    MM = matrix_space(FF, n, 1)
    return SymplecticVector{n, d}(MM(z), MM(x))
end


# Basic operations on symplectic vectors

"""
    hash(v::SymplecticVector{n, d}) -> UInt

Compute a hash for the given `SymplecticVector`.
"""
function hash(v::SymplecticVector{n, d}) where {n, d}
    return hash((v.z, v.x))
end

"""
    show(io::IO, v::SymplecticVector{n, d})

Print a clean representation of the `SymplecticVector` to the IO stream.
"""
function show(io::IO, v::SymplecticVector{n, d}) where {n, d}
    print(io, reduce(*, [
            "[",
            join(permutedims(hcat(string.(v.z[:,1]), repeat([" "], n)))[:])[1:end-1],
            "|",
            join(permutedims(hcat(string.(v.x[:,1]), repeat([" "], n)))[:])[1:end-1],
            "]"
        ]))
end

"""
    rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, d}}) -> SymplecticVector{n, d}

Generate a random symplectic vector using the given random number generator.
"""
function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, d}}) where {n, d}
    FF = finite_field(d)[1]
    return SymplecticVector{n, d}(rand(FF, n), rand(FF, n))
end

"""
    rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, d}}, dims...) -> Array{SymplecticVector{n, d}}

Generate an array of random symplectic vectors with the specified dimensions using the given random number generator.
"""
function rand(rng::AbstractRNG, ::SamplerType{SymplecticVector{n, d}}, dims...) where {n, d}
    FF = finite_field(d)[1]
    return SymplecticVector{n, d}.(rand(FF, n, dims...), rand(FF, n, dims...))
end


# Functions for extracting information from SymplecticVector types

"""
    dimension(v::SymplecticVector{n, d}) -> Int

Return the dimension (2n) of the symplectic vector.
"""
dimension(v::SymplecticVector{n, d}) where {n, d} = 2n

"""
    halfdimension(v::SymplecticVector{n, d}) -> Int

Return the half-dimension (n) of the symplectic vector.
"""
halfdimension(v::SymplecticVector{n, d}) where {n, d} = n


# Basic arithmetic with symplectic vectors

"""
    zero(::Type{SymplecticVector{n, d}}) -> SymplecticVector{n, d}

Return the zero symplectic vector.
"""
function zero(::Type{SymplecticVector{n, d}}) where {n, d}
    FF = finite_field(d)[1]
    return SymplecticVector{n, d}(repeat([zero(FF)], n), repeat([zero(FF)], n))
end

function zero(::SymplecticVector{n, d}) where {n, d}
    FF = finite_field(d)[1]
    return SymplecticVector{n, d}(repeat([zero(FF)], n), repeat([zero(FF)], n))
end


==(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = (u.z == v.z) && (u.x == v.x)

+(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = SymplecticVector{n, d}(u.z+v.z, u.x+v.x)
-(v::SymplecticVector{n, d}) where {n, d} = SymplecticVector{n, d}(-v.z, -v.x)
-(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = u + (-v)
*(k::FqFieldElem, v::SymplecticVector{n, d}) where {n, d} = SymplecticVector{n, d}(k * v.z, k * v.x)
*(v::SymplecticVector{n, d}, k::FqFieldElem) where {n, d} = SymplecticVector{n, d}(k * v.z, k * v.x)

"""
    dotproduct(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) -> FqFieldElem

Compute the dot product of two symplectic vectors.
"""
function dotproduct(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    (transpose(u.z) * v.z + transpose(u.x) * v.x)[1]
end
⋅(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = innerproduct(u, v)

"""
    symplecticform(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) -> FqFieldElem

Compute the symplectic form of two symplectic vectors.
"""
function symplecticform(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    (transpose(u.z) * v.x - transpose(u.x) * v.z)[1]
end
⋆ = symplecticform


###################################################################
##  Types and methods for subspaces of symplectic vector spaces  ##
###################################################################

"""
    islinearlyindependent(vectors::Vector{SymplecticVector{n, d}}) -> Bool

Check if a set of symplectic vectors is linearly independent.
"""
function islinearlyindependent(vectors::Vector{SymplecticVector{n, d}}) where {n, d}
    M = transpose(vcat(
        hcat(collect(vectors[i].z for i ∈ 1:length(vectors))...),
        hcat(collect(vectors[i].x for i ∈ 1:length(vectors))...)
        ))
    r, A = rref(M)
    if length(vectors) == r
        return true
    else
        return false
    end
end

"""
    isisotropic(vectors::Vector{SymplecticVector{n, d}}) -> Bool

Check if a set of symplectic vectors forms an isotropic subspace.
"""
function isisotropic(vectors::Vector{SymplecticVector{n, d}}) where {n, d}
    for i ∈ 1:length(vectors)
        for j ∈ (i+1):length(vectors)
            if !iszero(vectors[i] ⋆ vectors[j])
                return false
            end
        end
    end
    return true
end


abstract type AbstractSubspace end

"""
    Subspace{n, d}(basis::Vector{SymplecticVector{n, d}}; check::Bool=true)

Construct a subspace of a symplectic vector space with the provided basis.
Checks for linear independence of the basis vectors if `check=true`.
"""
struct Subspace{n, d} <: AbstractSubspace
    basis::Vector{SymplecticVector{n, d}}
    function Subspace{n, d}(basis::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
        if check
            if !islinearlyindependent(basis)
                error("Basis is not linearly independent")
            end
        end
        new(basis)
    end
end

function Subspace(basis::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
    return Subspace{n, d}(basis, check=check)
end

"""
    IsotropicSubspace{n, d}(basis::Vector{SymplecticVector{n, d}}; check::Bool=true)

Construct an isotropic subspace of a symplectic vector space with the provided basis.
Checks for isotropy and linear independence of the basis vectors if `check=true`.
"""
struct IsotropicSubspace{n, d} <: AbstractSubspace
    basis::Vector{SymplecticVector{n, d}}
    function IsotropicSubspace{n, d}(basis::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
        if check
            if !isisotropic(basis)
                error("Subspace is not isotropic")
            end
            if !islinearlyindependent(basis)
                error("Basis is not linearly independent")
            end
        end
        new(basis)
    end
end

function IsotropicSubspace(basis::Vector{SymplecticVector{n, d}}, check::Bool=true) where {n, d}
    return IsotropicSubspace{n, d}(basis, check=check)
end

"""
    LagrangianSubspace{n, d}(basis::Vector{SymplecticVector{n, d}}; check::Bool=true)

Construct a Lagrangian subspace of a symplectic vector space with the provided basis.
Checks for isotropy, linear independence, and that the dimension is `n` if `check=true`.
"""
struct LagrangianSubspace{n, d} <: AbstractSubspace
    basis::Vector{SymplecticVector{n, d}}
    function LagrangianSubspace{n, d}(basis::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
        if check
            if length(basis) ≠ n
                error("Not enough basis vectors")
            end
            if !isisotropic(basis)
                error("Subspace is not isotropic")
            end
            if !islinearlyindependent(basis)
                error("Basis is not linearly independent")
            end
        end
        new(basis)
    end
end

function LagrangianSubspace(basis::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
    return LagrangianSubspace{n, d}(basis, check=check)
end

function dimension(subspace::T) where {T<:AbstractSubspace}
    return length(subspace.basis)
end


#############################################
##  Types and methods for symplectic maps  ##
#############################################

function symplecticgrouporder(n, d)
    order = BigInt(d)^(n^2)
    for k ∈ 1:n
        order *= BigInt(d)^(2k) - 1
    end
    if order ≤ typemax(Int64)
        return Int64(order)
    elseif order ≤ typemax(Int128)
        return Int128(order)
    end
    return order
end

function issymplectic(z::Vector{SymplecticVector{n, d}}, x::Vector{SymplecticVector{n, d}}) where {n, d}
    for i ∈ 1:n
        for j ∈ (i+1):n
            if !iszero(z[i] ⋆ z[j]) || !iszero(x[i] ⋆ x[j])
                return false
            end
        end
    end
    for i ∈ 1:n
        for j ∈ 1:n
            if i == j
                if !isone(z[i] ⋆ x[j])
                    return false
                end
            else
                if !iszero(z[i] ⋆ x[j])
                    return false
                end
            end
        end
    end
    return true
end

abstract type AbstractSymplecticMap end

struct SymplecticMap{n, d}
    z_image::Vector{SymplecticVector{n, d}}
    x_image::Vector{SymplecticVector{n, d}}
    function SymplecticMap{n, d}(z_image::Vector{SymplecticVector{n, d}}, x_image::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
        if check
            if length(z_image) ≠ n || length(x_image) ≠ n
                error("Nonsquare matrix")
            end
            if !issymplectic(z_image, x_image)
                error("Map is not symplectic")
            end
        end

        new(z_image, x_image)
    end
end

function SymplecticMap(z_image::Vector{SymplecticVector{n, d}}, x_image::Vector{SymplecticVector{n, d}}; check::Bool=true) where {n, d}
    return SymplecticMap{n, d}(z_image, x_image; check=check)
end

function symplecticmap(map::SymplecticMap{n, d}, v::SymplecticVector{n, d}) where {n, d}
    matrix = hcat(
        vcat(
            hcat(collect(map.z_image[i].z for i ∈ 1:n)...),
            hcat(collect(map.z_image[i].x for i ∈ 1:n)...)
        ),
        vcat(
            hcat(collect(map.x_image[i].z for i ∈ 1:n)...),
            hcat(collect(map.x_image[i].x for i ∈ 1:n)...)
        )
    )
    data = Array(matrix * vcat(v.z, v.x))
    return SymplecticVector{n, d}(data[1:n], data[(n+1):2n])
end

*(map::SymplecticMap{n, d}, v::SymplecticVector{n, d}) where {n, d} = symplecticmap(map, v)


# Transvections
struct Transvection{n, d}
    h::SymplecticVector{n, d}
    λ::FqFieldElem
    function Transvection{n, d}(h::SymplecticVector{n, d}, λ::FqFieldElem) where {n, d}
        new(h, λ)
    end
end

function Transvection{n, d}(h::SymplecticVector{n, d}, λ::T) where {n, d, T<:Integer}
    FF = finite_field(d)[1]
    return Transvection{n, d}(h, FF(λ))
end

function Transvection(h::SymplecticVector{n, d}, λ) where {n, d}
    return Transvection{n, d}(h, λ)
end

function Transvection(h::SymplecticVector{n, d}) where {n, d}
    return Transvection{n, d}(h, 0)
end


function transvection(t::Transvection{n, d}, v::SymplecticVector{n, d}) where {n, d}
    return v + (t.λ * (v ⋆ t.h) * t.h)
end

function transvection(h::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    return transvection(Transvection(h, 0), v)
end

*(h::Transvection{n, d}, v::SymplecticVector{n, d}) where {n, d} = transvection(h, v)