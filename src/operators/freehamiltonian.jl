# free fermions



struct FreeFermionicHamiltonian <: AbstractFermionicHamiltonian
	data::Vector{TwoBodyTerm}
end

FreeFermionicHamiltonian() = FreeFermionicHamiltonian(Vector{TwoBodyTerm}())

Base.push!(x::FreeFermionicHamiltonian, m::TwoBodyTerm) = push!(x.data, m)
Base.copy(x::FreeFermionicHamiltonian) = FreeFermionicHamiltonian(copy(x.data))

change_positions(x::FreeFermionicHamiltonian, m::AbstractDict{Int, Int}) = FreeFermionicHamiltonian([change_positions(item, m) for item in x.data])


FermionicHamiltonian(x::FreeFermionicHamiltonian) = FermionicHamiltonian(convert(Vector{AbstractFermionicTerm}, x.data))

function FreeFermionicHamiltonian(x::FermionicHamiltonian)
	for item in x.data
		isa(item, TwoBodyTerm) || throw(ArgumentError("Can not convert an interacting hamiltonian into a free hamiltonian."))
	end
	return FreeFermionicHamiltonian(convert(Vector{TwoBodyTerm}, x.data))
end
free(x::FermionicHamiltonian) = FreeFermionicHamiltonian(x)


# """
# 	coefficient_matrix(x::FreeFermionicHamiltonian)

# Retuen the coefficient matrix for free fermion model
# """
# function coefficient_matrix(x::FreeFermionicHamiltonian, L::Int)
# 	@assert (length(x) <= L)
# 	@assert is_constant(x)
# 	m = spzeros(scalar_type(x), L, L)
# 	for item in x.data
# 		i, j = positions(item)
# 		m[i, j] = value(coeff(item))
# 	end
# 	return m
# end
# coefficient_matrix(x::FreeFermionicHamiltonian) = coefficient_matrix(x, length(x))


# """
# 	freefermion_stepper!(rhoout::AbstractMatrix, rho::AbstractMatrix, h::AbstractMatrix, c::Number)

# Return ρₒ += c (hρ - ρh)

# If h is the coefficient matrix of a free fermion Hamiltonian, one should take 
# the transpose of it as the input of this function for correct time evolution
# """
# function freefermion_stepper!(rhoout::AbstractMatrix, rho::AbstractMatrix, h::AbstractMatrix, c::Number)
# 	@assert size(rhoout) == size(rho) == size(h)
# 	mul!(rhoout, h, rho, c, one(c))
# 	mul!(rhoout, rho, h, -c, one(c))
# 	return rhoout
# end

# function freefermion_stepper(rho::AbstractMatrix, h::AbstractMatrix, c::Number)
# 	rhoout = zeros(promote_type(eltype(rho), eltype(h), typeof(c)), size(rho))
# 	return freefermion_stepper!(rhoout, rho, h, c)
# end