################################################
##  Types and methods for symplectic vectors  ##
################################################

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
⋅(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d} = innerproduct(u, v)

function symplecticform(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    (transpose(u.z) * v.x - transpose(u.x) * v.z)[1]
end
⋆ = symplecticform


###################################################################
##  Types and methods for subspaces of symplectic vector spaces  ##
###################################################################

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
    function SymplecticMap{n, d}(z_image, x_image; check::Bool=true) where {n, d}
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