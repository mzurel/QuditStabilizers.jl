############################################################################################
##  Tools for working with the stabilizer formalism on odd-prime-dimensional qudits.      ##
############################################################################################
module QuditStabilizers

using LinearAlgebra
using Kronecker
using Nemo

import Base: (==), (+), (-), (*), hash, isequal, iszero, show, rand, zero
import Random: AbstractRNG, SamplerType

export SymplecticVector, SymplecticSubspace
export dimension, halfdimension, data
export innerproduct, ⋅, symplecticform, ⋆

export Pauli
export compose, *, operator

include("Symplectic.jl")
include("Stabilizers.jl")

end
# End module