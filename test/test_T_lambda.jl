
@testset "λ T" begin
    for n_rule = 1:7
        # table = Radau.RadauTable(k)
        # @test table.T * diagm(0=>table.λ) * table.T⁻¹ ≈ inv(table.A)

        A, b, c, bi, b_hat = radau_butcher_table_plus(2 * n_rule - 1)
        inv_A = inv(A)
        λ, T = big_eigen(inv_A)

        @test T * diagm(0=>λ) * inv(T) ≈ inv_A

        if isodd(n_rule)
            @test λ[1] == real.(λ[1])
        end
    end
end
