#############################################
##  Types and methods for Pauli operators  ##
#############################################

struct Pauli{n, d}
    a::SymplecticVector{n, d}
    ϕ::FqFieldElem
    function Pauli{n, d}(a::SymplecticVector{n, d}, ϕ::FqFieldElem) where {n, d}
        new(a, ϕ)
    end
end

##  Outer constructors for convenience
function Pauli(a::SymplecticVector{n, d}, ϕ::FqFieldElem) where {n, d}
    return Pauli{n, d}(a, ϕ)
end

function Pauli(a::SymplecticVector{n, d}, ϕ::T) where {n, d, T<:Integer}
    FF = finite_field(d)[1]
    return Pauli{n, d}(a, FF(ϕ))
end

function Pauli(a::SymplecticVector{n, d}) where {n, d}
    return Pauli{n, d}(a, 0)
end

# Overloading builtin functions for SymplecticVector types
function hash(P::Pauli{n, d}) where {n, d}
    return hash((P.ϕ, P.a))
end

function show(io::IO, P::Pauli{n, d}) where {n, d}
    print(io, (P.a, P.ϕ))
end

# Functions for generating random Paulis
function rand(rng::AbstractRNG, ::SamplerType{Pauli{n, d}}) where {n, d}
    FF = finite_field(d)[1]
    return Pauli{n, d}(rand(SymplecticVector{n, d}), rand(FF))
end

function rand(rng::AbstractRNG, ::SamplerType{Pauli{n, d}}, dims...) where {n, d}
    FF = finite_field(d)[1]
    return Pauli{n, d}.(rand(SymplecticVector{n, d}, dims...), rand(FF, dims...))
end

## Basic operations on Pauli types
function compose(P::Pauli{n, d}, Q::Pauli{n, d}) where {n, d}
    return Pauli{n, d}(P.a + Q.a, P.ϕ + Q.ϕ + ((P.a ⋆ Q.a) / 2))
end
*(P::Pauli{n, d}, Q::Pauli{n, d}) where {n, d} = compose(P, Q)


function operator(P::Pauli{n, d}) where {n, d}
    ω = exp(2.0im * π / d)
    Z = diagm([ω^k for k ∈ 0:(d-1)])
    X = diagm(-1=>ones(ComplexF64, d-1), d-1=>ones(ComplexF64, 1))

    op = (Z ^ Integer(lift(ZZ, P.a.z[1]))) * (X ^ Integer(lift(ZZ, P.a.x[1])))
    for k ∈ 2:n
        op = op ⊗ ((Z ^ Integer(lift(ZZ, P.a.z[k]))) * (X ^ Integer(lift(ZZ, P.a.x[k]))))
    end
    op = ω^Integer(lift(ZZ, P.ϕ - ((transpose(P.a.z) * P.a.x)[1] / 2))) * op
    return op
end


struct Clifford{n, d}
    S::SymplecticMap{n, d}
    a::SymplecticVector{n, d}
    function Clifford{n, d}(S::SymplecticMap{n, d}, a::SymplecticVector{n, d}) where {n, d}
        new(S, a)
    end
end

function hash(C::Clifford{n, d}) where {n, d}
    return hash((C.S, C.a))
end

function rand(rng::AbstractRNG, ::SamplerType{Clifford{n, d}}) where {n, d}
    return Clifford{n, d}(rand(SymplecticMap{n, d}), rand(SymplecticVector{n, d}))
end

function *(C::Clifford{n, d}, P::Pauli{n, d}) where {n, d}
    return Pauli{n, d}(C.S * P.a, P.ϕ + (C.a ⋆ P.a))
end


struct StabilizerState{n, d}
    I::LagrangianSubspace{n, d}
    r::Vector{FqFieldElem}
end

function hash(σ::StabilizerState{n, d}) where {n, d}
    return hash((σ.I, σ.r))
end

function rand(rng::AbstractRNG, ::SamplerType{StabilizerState{n, d}}) where {n, d}
    return StabilizerState{n, d}(rand(LagrangianSubspace{n, d}), rand(finite_field(d)[1], n))
end
