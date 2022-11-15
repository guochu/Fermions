
using JSON

function read_data(pathname)
	data = JSON.parsefile(pathname)
	E0 = data["E0"]
	L = data["L"]
	t = data["t"]
	t = [t...]
	v = data["v"]
	v = [v...]
	return E0, reshape(t, (L, L)), reshape(v, (L, L, L, L))
end

read_lih_data() = read_data("lih.json")

const LiH_FCI_ENERGY = -7.78446028003123

function ed_lih_spinless()
	E0, t, v = read_lih_data()

	h = hamiltonian(SpinfulFermion(), t, v)

	mpo = MPO(h)
	# println(bond_dimensions(mpo))

	eigvalue, eigvector = ground_state(mpo, DMRG(D=40))

	energy = eigvalue + E0
	return abs(energy - LiH_FCI_ENERGY) < 1.0e-8
end

function ed_lih_spinful()
	E0, t, v = read_lih_data()

	h = hamiltonian(SpinlessFermion(), t, v)

	mpo = MPO(h)
	# println(bond_dimensions(mpo))

	eigvalue, eigvector = ground_state(mpo, DMRG(D=40))

	energy = eigvalue + E0
	return abs(energy - LiH_FCI_ENERGY) < 1.0e-8
end

println("check ground state of quantum chemistry system...")

@testset "check LiH ground state with ED" begin
	@test ed_lih_spinless()
	@test ed_lih_spinful() 
end

