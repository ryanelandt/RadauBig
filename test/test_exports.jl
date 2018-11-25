@testset "exports" begin
    # Ensure that every exported name is actually defined
    for name in names(RadauBig)
        @test isdefined(RadauBig, name)
        !isdefined(RadauBig, name) && println("MISSING: ", name)
    end
end
