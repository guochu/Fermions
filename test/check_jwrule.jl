

function check_spin_operators()
	cr = SpinCreation()
	an = SpinAnnihilation()
	n1 = SpinOccupationOne()
	n0 = SpinOccupationZero()
	jw = JWOperator()
	id = IdentityOperator()
	em = EmptyOperator()
	v1 = (cr * cr == em) && (cr * an == n1) && (cr * n1 == em) && (cr * n0 == cr) && (cr * jw == cr) && (cr * id == cr) && (cr * em == em)
	v2 = (an * cr == n0) && (an * an == em) && (an * n1 == an) && (an * n0 == em) && (an * jw == -an) && (an * id == an) && (an * em == em)
	v3 = (n1 * cr == cr) && (n1 * an == em) && (n1 * n1 == n1) && (n1 * n0 == em) && (n1 * jw == -n1) && (n1 * id == n1) && (n1 * em == em)
	v4 = (n0 * cr == em) && (n0 * an == an) && (n0 * n1 == em) && (n0 * n0 == n0) && (n0 * jw == n0) && (n0 * id == n0) && (n0 * em == em)
	v5 = (jw * cr == -cr) && (jw * an == an) && (jw * n1 == -n1) && (jw * n0 == n0) && (jw * jw == id) && (jw * id == jw) && (jw *em == em)
	v6 = (id * cr == cr) && (id * an == an) && (id * n1 == n1) && (id * n0 == n0) && (id * jw == jw) && (id * id == id) && (id * em == em)
	v7 = (em * cr == em) && (em * an == em) && (em * n1 == em) && (em * n0 == em) && (em * jw == em) && (em * id == em) && (em * em == em)
	v8 = (cr' == an) && (an' == cr) && (n1' == n1) && (n0' == n0) && (jw' == jw) && (id' == id) && (em' == em)
	return v1 && v2 && v3 && v4 && v5 && v6 && v7 && v8
end


function check_jw_rules()

	# single fermionic operator
	mf = SpinlessFermionicTerm([2], [FermionicCreation()])
	ms = jw_rule(mf)
	v1 = (ms.positions == [1,2]) && (ms.op[1] == JWOperator()) && (ms.op[2] == SpinCreation())

	mf = SpinlessFermionicTerm([3], [FermionicAnnihilation()])
	ms = jw_rule(mf)
	v2 = (ms.positions == [1,2,3]) && (ms.op[1] == JWOperator()) && (ms.op[2] == JWOperator()) && (ms.op[3] == SpinAnnihilation())

	mf = SpinlessFermionicTerm([5], [FermionicOccupation()])
	ms = jw_rule(mf)
	v3 = (ms.positions == [5]) && (ms.op[1] == SpinOccupationOne())

	# two fermionic operators
	mf = SpinlessFermionicTerm([2, 3], [FermionicCreation(), FermionicAnnihilation()])
	ms = jw_rule(mf)
	v4 = (ms.positions == [2,3]) && (ms.op[1] == SpinCreation()) && (ms.op[2] == SpinAnnihilation())


	mf = SpinlessFermionicTerm([2, 4], [FermionicAnnihilation(), FermionicCreation()])
	ms = jw_rule(mf)
	v5 = (ms.positions == [2,3, 4]) && (ms.op[1] == SpinAnnihilation()) && (ms.op[2] == JWOperator()) && (ms.op[3] == SpinCreation())


	mf = SpinlessFermionicTerm([2, 4], [FermionicOccupation(), FermionicOccupation()])
	ms = jw_rule(mf)
	v6 = (ms.positions == [2, 4]) && (ms.op[1] == SpinOccupationOne()) && (ms.op[2] == SpinOccupationOne())


	mf = SpinlessFermionicTerm([2, 2], [FermionicCreation(), FermionicAnnihilation()])
	ms = jw_rule(mf)
	v7 = (ms.positions == [2]) && (ms.op[1] == SpinOccupationOne()) 

	# # this behaviour may be changed
	# mf = SpinlessFermionicTerm([2, 2], [FermionicCreation(), FermionicCreation()])
	# ms = jw_rule(mf)
	# mf2 = SpinlessFermionicTerm([2, 2], [FermionicAnnihilation(), FermionicAnnihilation()])
	# ms2 = jw_rule(mf)
	# v8 = isnothing(ms) && isnothing(ms2)

	

	return v1 && v2 && v3 && v4 && v5 && v6 && v7 
end


@testset "jw rules" begin
	@test check_spin_operators()
	@test check_jw_rules() 
end