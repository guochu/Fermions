

abstract type AbstractFermionicTerm end
# abstract type AbstractFermionicTwoBodyTerm <: AbstractFermionicTerm end
# abstract type AbstractFermionicFourBodyTerm <: AbstractFermionicTerm end

QuantumSpins.positions(x::AbstractFermionicTerm) = x.positions
QuantumSpins.coeff(x::AbstractFermionicTerm) = x.coeff
QuantumSpins.is_constant(x::AbstractFermionicTerm) = is_constant(coeff(x))
QuantumSpins.scalar_type(x::AbstractFermionicTerm) = scalar_type(coeff(x))

"""
	struct TwoBodyTerm <: AbstractFermionicTerm

Fermionic twobody term
"""
struct TwoBodyTerm <: AbstractFermionicTerm
	positions::Tuple{Int, Int}
	coeff::AbstractCoefficient
end

TwoBodyTerm(pos::Tuple{Int, Int}; coeff::AllowedCoefficient=1.) = TwoBodyTerm(pos, Coefficient(coeff))
TwoBodyTerm(i::Int, j::Int; kwargs...) = TwoBodyTerm((i, j); kwargs...)

function Base.adjoint(x::TwoBodyTerm)
	i, j = positions(x)
	return TwoBodyTerm((j, i), coeff=conj(coeff(x)))
end

Base.copy(x::TwoBodyTerm) = TwoBodyTerm(copy(positions(x)), copy(coeff(x)))

Base.:*(s::TwoBodyTerm, m::AllowedCoefficient) = TwoBodyTerm(positions(s), coeff=coeff(s) * Coefficient(m))
Base.:*(m::AllowedCoefficient, s::AbstractFermionicTerm) = s * m
Base.:/(s::AbstractFermionicTerm, m::AllowedCoefficient) = s * (1 / Coefficient(m))
Base.:+(s::AbstractFermionicTerm) = s
Base.:-(s::AbstractFermionicTerm) = (-1) * s

"""
	struct FourBodyTerm <: AbstractFermionicTerm

Fermionic fourbody term
"""
struct FourBodyTerm <: AbstractFermionicTerm
	positions::NTuple{4, Int}
	coeff::AbstractCoefficient
end

FourBodyTerm(pos::NTuple{4, Int}; coeff::AllowedCoefficient=1.) = FourBodyTerm(pos, Coefficient(coeff))
FourBodyTerm(i::Int, j::Int, k::Int, l::Int; kwargs...) = FourBodyTerm((i, j, k, l); kwargs...)

Base.:*(s::FourBodyTerm, m::AllowedCoefficient) = FourBodyTerm(positions(s), coeff=coeff(s) * Coefficient(m))

function Base.adjoint(x::FourBodyTerm)
	i, j, k, l = positions(x)
	return fourbody((l,k,j,i), coeff=conj(coeff(x)))
end

function consolidate(x::TwoBodyTerm, symmetry::FermionicLatticeType) 
	i, j = positions(x)
	return twobody(i, j, coeff(x), symmetry)
end

function consolidate(x::FourBodyTerm, symmetry::FermionicLatticeType)
	i, j, k, l = positions(x)
	return fourbody(i,j,k,l, coeff(x), symmetry)
end


change_positions(x::TwoBodyTerm, m::AbstractDict{Int, Int}) = TwoBodyTerm(Tuple(m[k] for k in positions(x)), coeff=coeff(x))
change_positions(x::FourBodyTerm, m::AbstractDict{Int, Int}) = FourBodyTerm(Tuple(m[k] for k in positions(x)), coeff=coeff(x))



