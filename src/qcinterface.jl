
"""
	hamiltonian(lattice::FermionicLattice, t::Array{<:Real, 2}, v::Array{<:Real, 4}; atol::Real=1.0e-14)

Return the quantum chemistry Hamiltonian with the correct symmetry set in lattice, 
given the twobody coefficients t and fourbody coefficients v.

Definiton of the Hamiltonian can be found in docs/Hamiltonian.jpeg, for a molecule 
with L spatial orbital, t should have size (L, L) and v should have size (L,L,L,L).

Three different symmetries are supported for the Fermionic lattice:
* NoSymmetry, no symmetry
* ChargeCharge, U₁ ⊗ U₁ symmetry for up and down charge conservation
* SpinCharge, U₁ ⊗ SU₂ spin SU₂ symmetry and total U₁ chagre conservation
"""
function hamiltonian(lat::FermionicLatticeType, t::Array{<:Real, 2}, v::Array{<:Real, 4}; kwargs...)
	@assert size(t, 1) == size(t, 2) == size(v, 1) == size(v, 2) == size(v, 3) == size(v, 4)
	return consolidate(FermionicHamiltonian(t, v; kwargs...), lat)
end 


function FermionicHamiltonian(t::Array{<:Real, 2}, v::Array{<:Real, 4}; atol::Real=1.0e-14)
	@assert size(t, 1) == size(t, 2) == size(v, 1) == size(v, 2) == size(v, 3) == size(v, 4)
	L = size(t, 1)
	terms = AbstractFermionicTerm[]

	for pos1 in 1:L
		for pos2 in 1:L
			tmp = t[pos1, pos2]
			if abs(tmp) > atol
				push!(terms, TwoBodyTerm(pos1, pos2, coeff=tmp))
			end
		end
	end
	# println("number of twobody terms $(length(terms))")
	for pos1 in 1:L
		for pos2 in 1:L
			for pos3 in 1:L
				for pos4 in 1:L
					tmp = v[pos1, pos2, pos3, pos4]
					if abs(tmp) > atol
						push!(terms, FourBodyTerm(pos1, pos2, pos3, pos4, coeff=tmp))
					end
				end
			end
		end
	end
	isempty(terms) && error("vanishing hamiltonian")
	# println("total number of terms $(length(terms))")
	return FermionicHamiltonian(terms)
end

