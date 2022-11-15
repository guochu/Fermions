
abstract type AbstractSpinOperator end

struct SpinCreation <: AbstractSpinOperator end
struct SpinAnnihilation <: AbstractSpinOperator end
struct SpinOccupationOne <: AbstractSpinOperator end
struct SpinOccupationZero <: AbstractSpinOperator end
struct JWOperator <: AbstractSpinOperator end
struct IdentityOperator <: AbstractSpinOperator end
struct EmptyOperator <: AbstractSpinOperator end


const _sp = "+"
const _sm = "-"
const _one = "↑"
const _zero = "↓"
const _jw = "JW"
const _id = "I"
const _empty = "o"

const spin_operators = [_sp, _sm, _one, _zero, _jw, _id, _empty]

_repr(m::SpinCreation) = _sp
_repr(m::SpinAnnihilation) = _sm
_repr(m::SpinOccupationOne) = _one
_repr(m::SpinOccupationZero) = _zero
_repr(m::JWOperator) = _jw
_repr(m::IdentityOperator) = _id
_repr(m::EmptyOperator) = _empty

function SpinOperator(s::String)
	if s == _sp
		return SpinCreation()
	elseif s == _sm
		return SpinAnnihilation()
	elseif s == _one
		return SpinOccupationOne()
	elseif s == _zero
		return SpinOccupationZero()
	elseif s == _jw
		return JWOperator()
	elseif s == _id
		return IdentityOperator()
	elseif s == _empty
		return EmptyOperator()
	else
		throw(ArgumentError("fermionic operator must be in {+, -, ↑, ↓, JW, I, o}"))
	end
end


struct Negative{M <: AbstractSpinOperator}
	parent::M
end

Base.adjoint(m::SpinCreation) = SpinAnnihilation()
Base.adjoint(m::SpinAnnihilation)= SpinCreation()
Base.adjoint(m::Union{SpinOccupationOne, SpinOccupationZero, JWOperator, IdentityOperator, EmptyOperator}) = m
Base.adjoint(m::Negative) = Negative(m.parent')

Base.:-(x::Negative) = x.parent
Base.:-(x::AbstractSpinOperator) = Negative(x)
Base.:-(x::EmptyOperator) = x

Base.:*(x::AbstractSpinOperator, y::Negative) = Negative(x * y.parent)
Base.:*(x::Negative, y::AbstractSpinOperator) = Negative(x.parent * y)
Base.:*(x::Negative, y::Negative) = x.parent * y.parent



Base.:*(x::SpinCreation, y::SpinCreation) = EmptyOperator()
Base.:*(x::SpinCreation, y::SpinAnnihilation) = SpinOccupationOne()
Base.:*(x::SpinCreation, y::SpinOccupationOne) = EmptyOperator()
Base.:*(x::SpinCreation, y::SpinOccupationZero) = SpinCreation()
Base.:*(x::SpinCreation, y::JWOperator) = SpinCreation()



Base.:*(x::SpinAnnihilation, y::SpinCreation) = SpinOccupationZero()
Base.:*(x::SpinAnnihilation, y::SpinAnnihilation) = EmptyOperator()
Base.:*(x::SpinAnnihilation, y::SpinOccupationOne) = SpinAnnihilation()
Base.:*(x::SpinAnnihilation, y::SpinOccupationZero) = EmptyOperator()
Base.:*(x::SpinAnnihilation, y::JWOperator) = -x



Base.:*(x::SpinOccupationOne, y::SpinCreation) = SpinCreation()
Base.:*(x::SpinOccupationOne, y::SpinAnnihilation) = EmptyOperator()
Base.:*(x::SpinOccupationOne, y::SpinOccupationOne) = SpinOccupationOne()
Base.:*(x::SpinOccupationOne, y::SpinOccupationZero) = EmptyOperator()
Base.:*(x::SpinOccupationOne, y::JWOperator) = -x


Base.:*(x::SpinOccupationZero, y::SpinCreation) = EmptyOperator()
Base.:*(x::SpinOccupationZero, y::SpinAnnihilation) = SpinAnnihilation()
Base.:*(x::SpinOccupationZero, y::SpinOccupationOne) = EmptyOperator()
Base.:*(x::SpinOccupationZero, y::SpinOccupationZero) = SpinOccupationZero()
Base.:*(x::SpinOccupationZero, y::JWOperator) = SpinOccupationZero()



Base.:*(x::JWOperator, y::SpinCreation) = -y
Base.:*(x::JWOperator, y::SpinAnnihilation) = y
Base.:*(x::JWOperator, y::SpinOccupationOne) = y * x
Base.:*(x::JWOperator, y::SpinOccupationZero) = y * x
Base.:*(x::JWOperator, y::JWOperator) = IdentityOperator()


Base.:*(x::IdentityOperator, y::IdentityOperator) = y
Base.:*(x::IdentityOperator, y::AbstractSpinOperator) = y
Base.:*(x::AbstractSpinOperator, y::IdentityOperator) = x


Base.:*(x::EmptyOperator, y::EmptyOperator) = y
Base.:*(x::EmptyOperator, y::AbstractSpinOperator) = x
Base.:*(x::AbstractSpinOperator, y::EmptyOperator) = y

Base.:*(x::IdentityOperator, y::EmptyOperator) = y
Base.:*(x::EmptyOperator, y::IdentityOperator) = x

struct SpinTerm <: AbstractTerm
	positions::Vector{Int}
	op::Vector{AbstractSpinOperator}
	coeff::AbstractCoefficient
end
function SpinTerm(positions::Vector{Int}, op::Vector{<:AbstractSpinOperator}; coeff::AllowedCoefficient=1.)
	@assert positions == sort(positions)
	return SpinTerm(positions, convert(Vector{AbstractSpinOperator}, op), Coefficient(coeff))
end 


Base.adjoint(x::SpinTerm) = SpinTerm(positions(x), adjoint.(op(x)), coeff=conj(coeff(x)))

function Base.:*(x::SpinTerm, y::SpinTerm)
	coef = coeff(x) * coeff(y)
	pos, opx, opy = _coerce_spinterms(x, y)
	tmp = mult_spinterms(opx, opy)
	if isnothing(tmp)
		return nothing
	else
		r, _c = tmp
		coef *= _c
		pos, r = get_rid_of_identity(pos, r)
		return SpinTerm(pos, r, coeff=coef)
	end
end

Base.:*(x::SpinTerm, y::Union{Number, Coefficient}) = SpinTerm(positions(x), op(x), coeff=coeff(x)*y)
Base.:*(x::Union{Number, Coefficient}, y::SpinTerm) = y * x

function _coerce_spinterms(x::SpinTerm, y::SpinTerm)
    opx = op(x)
    opy = op(y)
    pos = positions(x)
    if !(positions(x) == positions(y))
    	new_pos = sort([Set(vcat(positions(x), positions(y)))...])
    	new_opx = []
    	new_opy = []
    	for pos in new_pos
    		pos_x = findfirst(a->a==pos, positions(x))
    		pos_y = findfirst(a->a==pos, positions(y))
    		if isnothing(pos_x) && !(isnothing(pos_y))
    			push!(new_opx, IdentityOperator())
    			push!(new_opy, opy[pos_y])
    		elseif !(isnothing(pos_x)) && isnothing(pos_y)
    			push!(new_opx, opx[pos_x])
    			push!(new_opy, IdentityOperator())
    		elseif !(isnothing(pos_x)) && !(isnothing(pos_y))
    			push!(new_opx, opx[pos_x])
    			push!(new_opy, opy[pos_y])
    		else
    			throw(ArgumentError("why here?"))
    		end
    	end
    	opx = [new_opx...]
    	opy = [new_opy...]
    	pos = new_pos
    end
    return pos, opx, opy
end


has_empty(r::Vector) = !isnothing(findfirst(x->x==EmptyOperator(), r))

function mult_spinterms(a::Vector{<:AbstractSpinOperator}, b::Vector{<:AbstractSpinOperator})
	@assert length(a) == length(b)
	r = [aj * bj for (aj, bj) in zip(a, b)]
	if has_empty(r)
		return nothing
	else
		coef = 1
		rc = AbstractSpinOperator[]
		for item in r
			if isa(item, Negative)
				coef = -coef
				push!(rc, item.parent)
			elseif isa(item, AbstractSpinOperator)
				push!(rc, item)
			else
				error("unexpected eltype $(typeof(item))")
			end
		end
		return rc, coef
	end
end

function get_rid_of_identity(pos::Vector{Int}, ops::Vector{<:AbstractSpinOperator})
	@assert !has_empty(ops)
	valid_pos = findall(x->x != IdentityOperator(), ops)
	@assert !isempty(valid_pos)
	return pos[valid_pos], ops[valid_pos]
end


