abstract type AbstractFermionicOperator end

struct FermionicCreation <: AbstractFermionicOperator end
struct FermionicAnnihilation <: AbstractFermionicOperator end
struct FermionicOccupation <: AbstractFermionicOperator end


Base.adjoint(m::FermionicCreation) = FermionicAnnihilation()
Base.adjoint(m::FermionicAnnihilation) = FermionicCreation()
Base.adjoint(m::FermionicOccupation) = FermionicOccupation()


_repr(m::FermionicCreation) = "c†"
_repr(m::FermionicAnnihilation) = "c"
_repr(m::FermionicOccupation) = "c†c"

Base.show(io::IO, a::AbstractFermionicOperator) = print(io, "FermionicOperator(", _repr(a) ,")")

function FermionicOperator(s::String)
	if s == "c†"
		return FermionicCreation()
	elseif s == "c"
		return FermionicAnnihilation()
	elseif s == "c†c"
		return FermionicOccupation()
	else
		throw(ArgumentError("fermionic operator must be in {c†, c, c†c}"))
	end
end

QuantumSpins.positions(x::AbstractTerm) = x.positions
QuantumSpins.op(x::AbstractTerm) = x.op
QuantumSpins.coeff(x::AbstractTerm) = x.coeff
QuantumSpins.is_constant(x::AbstractTerm) = is_constant(coeff(x))


const AllowedCoefficient = Union{Number, Function, Coefficient}


struct SpinlessFermionicTerm <: AbstractTerm
	positions::Vector{Int}
	op::Vector{AbstractFermionicOperator}
	coeff::AbstractCoefficient
end

SpinlessFermionicTerm(positions::Vector{Int}, op::Vector{<:AbstractFermionicOperator}; coeff::AllowedCoefficient=1.) = SpinlessFermionicTerm(
	positions, convert(Vector{AbstractFermionicOperator}, op), Coefficient(coeff))

Base.adjoint(x::SpinlessFermionicTerm) = SpinlessFermionicTerm(positions(x), adjoint.(op(x)), coeff=conj(coeff(x)))