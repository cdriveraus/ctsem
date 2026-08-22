# These are tiny exported smoke-test helpers used by external callers, so the
# tests only check their public return behavior.
@testset "Core API" begin
    @test ContinuousTimeSEM.scalar_square(-3) == 9
    @test ContinuousTimeSEM.scalar_square(1.5) ≈ 2.25 atol=1e-12 rtol=1e-12
end
