############################################################################################
##  Tools for working with the stabilizer formalism on odd-prime-dimensional qudits.      ##
############################################################################################
module QuditStabilizers

using Base.Iterators
using LinearAlgebra
using Kronecker
using Nemo

import Base: (==), (+), (-), (*), hash, isequal, iszero, show, rand, zero
import Random: AbstractRNG, SamplerType

export SymplecticVector, SymplecticSubspace
export dimension, halfdimension, data, extend, extendfront
export innerproduct, (⋅), symplecticform, (⋆)

export Subspace, IsotropicSubspace, LagrangianSubspace
export islinearlyindependent, isisotropic

export SymplecticMap, Transvection
export symplecticgrouporder, issymplectic, symplecticmap, (*)
export transvection, findtransvection


export Pauli
export compose, operator

include("Symplectic.jl")
include("Stabilizers.jl")

end
# End module