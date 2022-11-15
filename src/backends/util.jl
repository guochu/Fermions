
function _dense_spin_matrices()
	p = spin_half_matrices()
	sp, sm, z = p["+"], p["-"], p["z"]
	I2 = one(z)
	n = (z + I2) / 2
	JW = -z
	return Dict(_sp=>sp, _sm=>sm, _id=>I2, _one=>n, _zero=>sm * sp, _jw=>JW, _empty=>zero(z))
end

const p_nosymmetry = _dense_spin_matrices()


"""
	abelian_term(m::SpinlessFermionicTerm, spins::Vector{<:SpinHalf}, lat::SpinlessFermion)

Return QTerm with spinless fermionic site (d=2)
"""
function abelian_term(m::SpinlessFermionicTerm, spins::Vector{<:SpinHalf}, lat::SpinlessFermion)
	mf = add_spin(m, spins)
	ms = jw_rule(mf)
	isnothing(ms) && return nothing
	p = p_nosymmetry
	ops = [p[_repr(x)] for x in op(ms)]
	return QTerm(positions(ms), ops, coeff=coeff(ms))
end


"""
	abelian_term(m::SpinlessFermionicTerm, spins::Vector{<:SpinHalf}, lat::SpinlessFermion)

Return QTerm with spinful fermionic site (d=4)
"""
function abelian_term(m::SpinlessFermionicTerm, spins::Vector{<:SpinHalf}, lat::SpinfulFermion)
	mf = add_spin(m, spins)
	ms = jw_rule(mf)
	isnothing(ms) && return nothing
	p = p_nosymmetry
	ops = [p[_repr(x)] for x in op(ms)]
	pos, ops = _to_spinful_ops(positions(ms), ops, p)
	return QTerm(pos, ops, coeff=coeff(ms))
end

function _to_spinful_ops(pos::Vector{Int}, v::Vector, p)
	@assert length(pos) == length(v)
	pos_start = div(pos[1]+1, 2)
	pos_end = div(pos[end]+1, 2)
	new_pos = Int[]
	new_v = []
	for i in pos_start:pos_end
		up = 2*i - 1
		down = 2 * i
		up_pos = findfirst(x -> x == up, pos)
		up_op = isnothing(up_pos) ? p[_id] : v[up_pos]
		down_pos = findfirst(x -> x == down, pos)
		down_op = isnothing(down_pos) ? p[_id] : v[down_pos]
		if !(isnothing(up_pos) && isnothing(down_pos))
			push!(new_pos, i)
			push!(new_v, kron(up_op, down_op))
		end
	end
	return new_pos, [new_v...]
end

add_spin(m::SpinlessFermionicTerm, spins::Vector{<:SpinHalf}) = SpinlessFermionicTerm(_add_spin(positions(m), spins), op(m), coeff=coeff(m))

function _add_spin(pos::Vector{Int}, spins::Vector{<:SpinHalf})
	@assert length(pos) == length(spins)
	new_pos = [2*a - 1 + _spin_pos(b) for (a, b) in zip(pos, spins)]
	return new_pos
end

_spin_pos(s::SpinUp) = 0
_spin_pos(s::SpinDown) = 1
