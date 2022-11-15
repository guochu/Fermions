module Fermions

using Reexport
using QuantumSpins
@reexport using QuantumSpins
import QuantumSpins


# fermionic site type
export FermionicLatticeType, SpinfulFermion, SpinlessFermion

# spinless definition and jordan-wigner transformation
export SpinHalf, SpinUp, SpinDown
export AbstractFermionicOperator, FermionicCreation, FermionicAnnihilation, FermionicOccupation, FermionicOperator, SpinlessFermionicTerm
export AbstractSpinOperator, SpinCreation, SpinAnnihilation, SpinOccupationOne, SpinOccupationZero, JWOperator, IdentityOperator, EmptyOperator
export SpinOperator, SpinTerm, Negative, jw_rule

# backends
export twobody, fourbody, creation

# fermionic hamiltonian wrapper
export AbstractFermionicTerm, TwoBodyTerm, FourBodyTerm, consolidate, change_positions, AbstractFermionicHamiltonian, FermionicHamiltonian, hamiltonian
export FreeFermionicHamiltonian, free, coefficient_matrix


abstract type AbstractTerm end
abstract type AbstractFermionicLattice end


# definition of fermionic site type
include("latticetype.jl")

# generating Hamiltonian without explitly preserving the non-Abelian symmetry
include("spinless/spin.jl")
include("spinless/jw/fermionicoperators.jl")
include("spinless/jw/spinoperators.jl")
include("spinless/jw/jwrule.jl")

# backends for genenrating twobody and fourbody terms
include("backends/backends.jl")


# hamiltonian wrapper
include("operators/fterm.jl")
include("operators/abstractdefs.jl")
include("operators/hamiltonian.jl")
include("operators/freehamiltonian.jl")

# interface for quantum chemistry
include("qcinterface.jl")


end

