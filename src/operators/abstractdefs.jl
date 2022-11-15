abstract type AbstractFermionicHamiltonian end

QuantumSpins.qterms(x::AbstractFermionicHamiltonian) = x.data
Base.isempty(x::AbstractFermionicHamiltonian) = isempty(x.data)
Base.length(x::AbstractFermionicHamiltonian) = maximum_site(x)

function QuantumSpins.scalar_type(x::AbstractFermionicHamiltonian)
	T = Float64
	for item in x.data
		T = promote_type(T, scalar_type(item))
	end
	return T
end
function QuantumSpins.is_constant(x::AbstractFermionicHamiltonian)
	for item in x.data
		is_constant(item) || return false
	end
	return true
end

function maximum_site(x::AbstractFermionicHamiltonian)
	isempty(x) && error("empty hamiltonian")
	m = 1
	for item in x.data
		m = max(m, maximum(positions(item)))
	end
	return m
end

function consolidate(x::AbstractFermionicHamiltonian, symmetry::FermionicLatticeType)
	terms = []
	for item in x.data
		_add!(terms, consolidate(item, symmetry))
	end
	isempty(terms) && error("vanishing hamiltonian")
	return QuantumOperator([terms...])
end


_add!(terms::Vector, m::QTerm) = push!(terms, m)
_add!(terms::Vector, m::Vector) = append!(terms, m)

