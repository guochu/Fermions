include("util.jl")

function twobody(pos1::Int, pos2::Int, c::AllowedCoefficient, lat::FermionicLatticeType)
	f = SpinlessFermionicTerm([pos1, pos2], [FermionicCreation(), FermionicAnnihilation()], coeff=c)
	r = []
	t = abelian_term(f, [SpinUp(), SpinUp()], lat)
	!isnothing(t) && push!(r, t)
	t = abelian_term(f, [SpinDown(), SpinDown()], lat)
	!isnothing(t) && push!(r, t)
	return r
end

function fourbody(pos1::Int, pos2::Int, pos3::Int, pos4::Int, v::AllowedCoefficient, lat::FermionicLatticeType)
	pos = [pos1, pos2, pos3, pos4]
	f = SpinlessFermionicTerm(pos, [FermionicCreation(), FermionicCreation(), FermionicAnnihilation(), FermionicAnnihilation()], coeff=0.5*Coefficient(v))
	r = []
	t = abelian_term(f, [SpinUp(), SpinUp(), SpinUp(), SpinUp()], lat)
	!isnothing(t) && push!(r, t)
	t = abelian_term(f, [SpinUp(), SpinDown(), SpinDown(), SpinUp()], lat)
	!isnothing(t) && push!(r, t)
	t = abelian_term(f, [SpinDown(), SpinUp(), SpinUp(), SpinDown()], lat)
	!isnothing(t) && push!(r, t)
	t = abelian_term(f, [SpinDown(), SpinDown(), SpinDown(), SpinDown()], lat)
	!isnothing(t) && push!(r, t)
	return r
end

function creation(pos::Int, spin::SpinHalf, lat::FermionicLatticeType)
	f = SpinlessFermionicTerm([pos], [FermionicCreation()])
	return abelian_term(f, [spin], lat)
end 





