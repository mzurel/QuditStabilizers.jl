#################################
##  Types for Pauli operators  ##
#################################

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


## Basic operations on Pauli types
function compose(P::Pauli{n, d}, Q::Pauli{n, d}) where {n, d}
    return Pauli{n, d}(P.a + Q.a, P.ϕ + Q.ϕ + ((P.a ⋆ Q.a) / 2))
end
*(P::Pauli{n, d}, Q::Pauli{n, d}) where {n, d} = compose(P, Q)


function operator(P::Pauli{n, d}) where {n, d}
    ω = exp(2.0im * π / 3)
    Z = diagm([ω^k for k ∈ 0:(d-1)])
    X = diagm(-1=>ones(ComplexF64, d-1), d-1=>ones(ComplexF64, 1))

    op = (Z ^ Integer(lift(ZZ, P.a.z[1]))) * (X ^ Integer(lift(ZZ, P.a.x[1])))
    for k ∈ 2:n
        op = op ⊗ ((Z ^ Integer(lift(ZZ, P.a.z[k]))) * (X ^ Integer(lift(ZZ, P.a.x[k]))))
    end
    op = ω^Integer(lift(ZZ, P.ϕ - ((transpose(P.a.z) * P.a.x)[1] / 2))) * op
    return op
end