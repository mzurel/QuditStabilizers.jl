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

"""
    extend(v::SymplecticVector{n, d}, m::Integer) where {n, d} -> SymplecticVector{m + n, d}

Embed a vector ``v∈ℤₚ²ⁿ`` in ``ℤₚ²⁽ᵐ⁺ⁿ⁾``.
"""
function extend(v::SymplecticVector{n, d}, m::Integer) where {n, d}
    MM = matrix_space(finite_field(d)[1], m, 1)
    return SymplecticVector{m + n, d}(vcat(v.z, zero(MM)), vcat(v.x, zero(MM)))
end

function extend(v::SymplecticVector{n, d}) where {n, d}
    return extend(v, 1)
end

function extendfront(v::SymplecticVector{n, d}, m::Integer) where{n, d}
    MM = matrix_space(finite_field(d)[1], m, 1)
    return SymplecticVector{m + n, d}(vcat(zero(MM), v.z), vcat(zero(MM), v.x))
end

function extendfront(v::SymplecticVector{n, d}) where {n, d}
    return extendfront(v, 1)
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

"""
    SymplecticMap{n, d}(z_image::Vector{SymplecticVector{n, d}}, x_image::Vector{SymplecticVector{n, d}}; check::Bool=true)

Construct a symplectic map from images of the symplectic basis vectors.
Checks if the map is symplectic if `check=true`.

# Arguments
    z_image::Vector{SymplecticVector{n, d}}: Vector of images for the z components.
    x_image::Vector{SymplecticVector{n, d}}: Vector of images for the x components.
    check::Bool=true: Whether to check if the map is symplectic.
"""
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

"""
    symplecticmap(map::SymplecticMap{n, d}, v::SymplecticVector{n, d})

Apply a SymplecticMap to a SymplecticVector, returning the image of the vector.

# Arguments
    map::SymplecticMap{n, d}: The symplectic map.
    v::SymplecticVector{n, d}: The vector to be mapped.
"""
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

"""
    symplecticmap(maps::Vector, v)

If a vector of SymplecticMaps `map` is given, the maps are applied to `v` sequentially.
"""
function symplecticmap(map::Vector, v)
    o = v
    for m ∈ maps
        o = symplecticmap(m, o)
    end
    return o
end

######################################################
##  Types and methods for symplectic transvections  ##
######################################################

"""
    Transvection{n, d}(h::SymplecticVector{n, d}, λ::FqFieldElem)

Construct a `Transvection` from a symplectic vector `h` and a scalar `λ`.

# Arguments
- `h::SymplecticVector{n, d}`: The symplectic vector that defines the transvection.
- `λ::FqFieldElem`: The scalar multiplier for the transvection.
"""
struct Transvection{n, d}
    h::SymplecticVector{n, d}
    λ::FqFieldElem
    function Transvection{n, d}(h::SymplecticVector{n, d}, λ::FqFieldElem) where {n, d}
        new(h, λ)
    end
end

""" Transvection{n, d}(h::SymplecticVector{n, d}, λ::T) where {T<Integer}

Outer constructor for Transvection that accepts an integer λ and converts it to a finite field element.

# Arguments
 - h::SymplecticVector{n, d}: The symplectic vector that defines the transvection.
 - λ::T<:Integer: An integer that will be converted to a finite field element.
"""
function Transvection{n, d}(h::SymplecticVector{n, d}, λ::T) where {n, d, T<:Integer}
    FF = finite_field(d)[1]
    return Transvection{n, d}(h, FF(λ))
end

"""
    Transvection(h::SymplecticVector{n, d}, λ) where {n, d}

Construct a Transvection from a symplectic vector h and a scalar λ. This function is a convenience constructor that infers the type parameters.

# Arguments
 - h::SymplecticVector{n, d}: The symplectic vector that defines the transvection.
 - λ: The scalar multiplier for the transvection.
"""
function Transvection(h::SymplecticVector{n, d}, λ) where {n, d}
    return Transvection{n, d}(h, λ)
end

"""
    Transvection(h::SymplecticVector{n, d}) where {n, d}

If no scalar multipler is given, assume it is 1.
"""
function Transvection(h::SymplecticVector{n, d}) where {n, d}
    return Transvection{n, d}(h, 1)
end

"""
    transvection(t::Transvection{n, d}, v::SymplecticVector{n, d}) where {n, d}

Apply a transvection t to a symplectic vector v.

# Arguments
    t::Transvection{n, d}: The transvection to be applied.
    v::SymplecticVector{n, d}: The symplectic vector to which the transvection is applied.
"""
function transvection(t::Transvection{n, d}, v::SymplecticVector{n, d}) where {n, d}
    return v + (t.λ * (v ⋆ t.h) * t.h)
end

function transvection(h::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    return transvection(Transvection(h, 1), v)
end

*(h::Transvection{n, d}, v::SymplecticVector{n, d}) where {n, d} = transvection(h, v)

"""
    transvection(t::Vector, v)

If a vector of transvections `t` is given, the transvections are applied to `v` sequentially.
"""
function transvection(ts::Vector, v)
    o = v
    for t ∈ ts
        o = transvection(t, o)
    end
    return o
end

"""
    findtransvection(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}

Given symplectic vectors `u` and `v`, finds a transvection `t` or a pair of transvections
`[t₁,t₂]` such that these transvections applied to `u` give `v`.
"""
function findtransvection(u::SymplecticVector{n, d}, v::SymplecticVector{n, d}) where {n, d}
    if u == v
        return Transvection(zero(SymplecticVector{n, d}))
    elseif !iszero(u ⋆ v)
        return Transvection(v - u, inv(u ⋆ v))
    end

    FF = finite_field(d)[1]
    for i ∈ 1:n
        if (!iszero(u.z[i]) || !iszero(u.x[i])) && (!iszero(v.z[i]) || !iszero(v.x[i]))
            for (a,b) ∈ product(repeat([[FF(0), FF(1)]], 2)...)
                z = zeros(FF, n); x = zeros(FF, n)
                z[i] = a; x[i] = b
                w = SymplecticVector{n, d}(z, x)
                if !iszero(u ⋆ w) && !iszero(v ⋆ w)
                    return [Transvection{n, d}(w - u, inv(u ⋆ w)), Transvection{n, d}(v - w, inv(w ⋆ v))]
                end
            end
        end
    end

    for i ∈ 1:n
        if !iszero(u.z[i]) || !iszero(u.x[i])
            for j ∈ 1:n
                if !iszero(v.z[j]) || !iszero(v.x[j])
                    for (a₁,b₁,a₂,b₂) ∈ product(repeat([[FF(0), FF(1)]], 4)...)
                        z = zeros(FF, n); x = zeros(FF, n)
                        z[i] = a₁; x[i] = b₁
                        z[j] = a₂; x[j] = b₂
                        w = SymplecticVector{n, d}(z, x)
                        if !iszero(u ⋆ w) && !iszero(v ⋆ w)
                            return [Transvection{n, d}(w - u, inv(u ⋆ w)), Transvection{n, d}(v - w, inv(w ⋆ v))]
                        end
                    end
                end
            end
        end
    end
end

function SYMPLECTICImproved(n, d, i)
    FF = finite_field(d)[1]
    MM = matrix_space(FF, n, 1)

    # Step 1
    s = BigInt(d)^(2n) - 1
    k = digits((i % s) + 1, base = d, pad = 2n)

    # Step 2
    e₁ = SymplecticVector{n, d}(zero(MM) + 1, zero(MM))
    f₁ = SymplecticVector{n, d}(k[1:n], k[(n+1):2n])

    # Step 3
    T = findtransvection(e₁, f₁)

    # Step 4
    b = FF(i ÷ s)
    b⃗ = digits(i ÷ s, base = d)[2:(2n-1)]

    # Step 5
    ẽ₁ = zero(MM) + 1
    ẽ₂ = zero(MM)
    for j ∈ 2:n
        ẽ₁[j] = b⃗[2j-3]
        ẽ₂[j] = b⃗[2j-2]
    end
    h₀ = transvection(T, SymplecticVector{n, d}(ẽ₁, ẽ₂))

    # Step 6
    T̃ = [Transvection(h₀), Transvection(b * f₁)]

    e₂ = SymplecticVector{n, d}(zero(MM), zero(MM) + 1)
    f₂ = transvection(vcat(T, T̃), e₂)

    # Step 7
    if n == 1
        return SymplecticMap{n, d}([f₁], [f₂])
    else
        B₍ₙ₋₁₎ = SYMPLECTICImproved(n-1, d, BigInt(floor((i ÷ s) / BigInt(d)^(2n-1))))
        Bₙ = SymplecticMap{n, d}(
            [transvection(vcat(T, T̃), v) for v ∈ [e₁, extendfront.(B₍ₙ₋₁₎.z_image)...]],
            [transvection(vcat(T, T̃), v) for v ∈ [e₂, extendfront.(B₍ₙ₋₁₎.x_image)...]]
        )
        return Bₙ
    end
end