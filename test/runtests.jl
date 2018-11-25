using Test
using RadauBig
using LinearAlgebra
using Polynomials: Poly, polyval, polyder, coeffs
using PolynomialRoots: roots


include("test_exports.jl")
include("test_T_lambda.jl")
include("test_generate_butcher_table.jl")


# TODO: add tests for calculating of Eigencalculatino of T and λ
