




struct FermionicHamiltonian <: AbstractFermionicHamiltonian
	data::Vector{AbstractFermionicTerm}

function FermionicHamiltonian(data::Vector{<:AbstractFermionicTerm})
	new(convert(Vector{AbstractFermionicTerm}, data))
end

end

FermionicHamiltonian() = FermionicHamiltonian(Vector{AbstractFermionicTerm}())

Base.push!(x::FermionicHamiltonian, m::AbstractFermionicTerm) = push!(x.data, m)
Base.copy(x::FermionicHamiltonian) = FermionicHamiltonian(copy(x.data))


change_positions(x::FermionicHamiltonian, m::AbstractDict{Int, Int}) = FermionicHamiltonian([change_positions(item, m) for item in x.data])


Base.:*(s::FermionicHamiltonian, m::AllowedCoefficient) = FermionicHamiltonian([item * m for item in s.data])
Base.:*(m::AllowedCoefficient, s::AbstractFermionicHamiltonian) = s * m
Base.:/(s::AbstractFermionicHamiltonian, m::AllowedCoefficient) = s * (1/m)
Base.:+(s::AbstractFermionicHamiltonian) = s
Base.:-(s::AbstractFermionicHamiltonian) = (-1) * s


Base.:+(a::AbstractFermionicTerm, b::AbstractFermionicTerm) = FermionicHamiltonian([a, b])
Base.:+(a::FermionicHamiltonian, b::FermionicHamiltonian) = FermionicHamiltonian(vcat(a.data, b.data))
function Base.:+(a::FermionicHamiltonian, b::AbstractFermionicTerm)
	ac = copy(a.data)
	push!(ac, b)
	return FermionicHamiltonian(ac)
end
Base.:+(a::AbstractFermionicTerm, b::FermionicHamiltonian) = b + a

