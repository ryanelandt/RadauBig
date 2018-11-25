
function save_radau_table_to_file(n_rule::Int64)
    @eval begin  # https://stackoverflow.com/questions/26071317/declaring-top-level-variables-in-julia-using-metaprogramming/26071597#26071597
        n_rule = $n_rule
        A, b, c, bi, b_hat = radau_butcher_table_plus(2 * n_rule - 1)
        inv_A = inv(A)
        lambda, T = big_eigen(inv_A)
        inv_T = inv(T)
        f_64_type(v::VecOrMat{BigFloat}) = Float64
        f_64_type(v::VecOrMat{Complex{BigFloat}}) = Complex{Float64}
        symbols_save = (:A, :inv_A, :b, :c, :lambda, :T, :inv_T, :bi, :b_hat)
        symbols_type = f_64_type.(eval.(symbols_save))

        ### Table Data
        radau_table_path = get_radau_table_path()
        save_path = joinpath(radau_table_path, string(n_rule) * "_rule")
        for (k, var) = enumerate(symbols_save)
            file_path = joinpath(save_path, string(var) * ".txt")
            open(file_path, "w") do io
                d_type_k = symbols_type[k]
                writedlm(io, d_type_k.(eval(var)))
            end
        end
    end

    ### Table Types
    save_path_types = joinpath(radau_table_path, "type.txt")
    open(save_path_types, "w") do io
        writedlm(io, string.(symbols_type))
    end

    ### Table Symbols
    save_path_names = joinpath(radau_table_path, "symbols.txt")
    open(save_path_names, "w") do io
        writedlm(io, string.(symbols_save))
    end
end
