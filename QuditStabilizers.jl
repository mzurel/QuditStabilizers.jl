############################################################################################
##  Tools for working with the stabilizer formalism on odd-prime-dimensional qudits.      ##
############################################################################################
module QuditStabilizers

using LinearAlgebra
using Kronecker
using Nemo

import Base: (==), (+), (-), (*), hash, isequal, iszero, show, rand, zero

export SymplecticVector, SymplecticSubspace
export dimension, halfdimension, data
export innerproduct, ⋅, symplecticform, ⋆

include("Symplectic.jl")

end
# End module