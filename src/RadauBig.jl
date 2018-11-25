module RadauBig

using ForwardDiff
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using Polynomials: Poly, polyval, polyder, coeffs
using PolynomialRoots: roots
using GenericLinearAlgebra: eigvals
using GenericSVD: svd

include("big.jl")
include("generate_butcher_table.jl")
include("save_radau_tables.jl")

export
    # big.jl
    big_eigen,

    # generate_butcher_table.jl
    find_real_eigenvalue,
    radau_butcher_table_plus,

    # save_radau_tables.jl
    save_radau_table_to_file

end
