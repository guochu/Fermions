



function jw_rule(m::SpinlessFermionicTerm) 
	@assert is_constant(m)
	ms = [_jw_rule_single(a, b) for (a, b) in zip(positions(m), op(m))]
	@assert !isempty(ms)
	if length(ms) == 1
		return ms[1] * value(coeff(m))
	else
		r = ms[1] * ms[2]
		isnothing(r) && return r
		for i in 3:length(ms)
			r = r * ms[i]
			isnothing(r) && return r
		end
		return r * value(coeff(m))
	end
end


function _jw_rule_single(pos::Int, op::FermionicCreation)
	_pos = collect(1:pos)
	_ops = Vector{AbstractSpinOperator}(undef, length(_pos))
	for i in 1:length(_pos)-1
		_ops[i] = JWOperator()
	end
	_ops[end] = SpinCreation()
	return SpinTerm(_pos, _ops)
end

function _jw_rule_single(pos::Int, op::FermionicAnnihilation)
	_pos = collect(1:pos)
	_ops = Vector{AbstractSpinOperator}(undef, length(_pos))
	for i in 1:length(_pos)-1
		_ops[i] = JWOperator()
	end
	_ops[end] = SpinAnnihilation()
	return SpinTerm(_pos, _ops)
end

_jw_rule_single(pos::Int, op::FermionicOccupation) = SpinTerm([pos], [SpinOccupationOne()])



