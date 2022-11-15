push!(LOAD_PATH, "/Users/guochu/Documents/QuantumSimulator/QuantumSpins/src")
push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")

using Test
using Fermions


include("check_jwrule.jl")
include("check_qc.jl")
