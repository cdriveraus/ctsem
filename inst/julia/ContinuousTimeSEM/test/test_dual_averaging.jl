# The dual-averaging restart, `_dual_restart!` (sample_adapt.jl), against the
# historical bug its own comment describes.
#
# `_dual_restart!` is called by the divergence-triggered `target_accept`
# auto-raise in `_run_chain` (sample_run.jl:291-319), which can fire on the
# very last warmup iteration -- and a restart there is followed by no further
# `_dual_update!` call before sampling reads the averaged step size back out.
# Before the fix documented at sample_adapt.jl:73-96, `logepsbar` was left at
# whatever an unstarted `_DualAverage` initialises it to, which for a restart
# with no prior updates at all is `log(0) -> eps = exp(0) = 1`; one chain in a
# real run took R-hat to 4.11 sampling at step size 1 against an adapted step
# size around 0.03. Grep across this suite finds no test that calls
# `_dual_restart!` or drives a real divergence-rate raise, so the exact path
# that produced that failure has no regression test.
#
# This pins the fix at the unit level, per the mechanism's own docstring:
# `_dual_restart!` sets `logepsbar = log(eps)` rather than leaving it at
# whatever the averaging held, so `_dual_final` after a restart with no
# further update returns the restarted step size, not 1.

@testset "a restart with no further update keeps the restarted step size, not 1" begin
    da = ContinuousTimeSEM._DualAverage(0.05, 0.8)
    # Warm it up first, as a real chain would, so the restart is genuinely
    # overriding accumulated averaging state rather than starting from
    # nothing.
    for accept in (0.6, 0.55, 0.7, 0.4)
        ContinuousTimeSEM._dual_update!(da, accept)
    end
    @test ContinuousTimeSEM._dual_final(da) != 1.0

    restarted_eps = 0.031
    ContinuousTimeSEM._dual_restart!(da, restarted_eps)
    # The historical bug, pinned directly: no `_dual_update!` call happens
    # between the restart and reading the step size back out, which is
    # exactly the "raise on the last warmup iteration" case.
    @test ContinuousTimeSEM._dual_final(da) ≈ restarted_eps
    @test ContinuousTimeSEM._dual_final(da) != 1.0
end

@testset "a restart from a fresh, never-updated average also avoids 1" begin
    # The other way the bug was reachable: `nwarmup = 0`, so no `_dual_update!`
    # call ever happens at all and `_init_stepsize`'s answer -- what a restart
    # would also centre on -- is what sampling must use.
    da = ContinuousTimeSEM._DualAverage(0.05, 0.8)
    restarted_eps = 0.4
    ContinuousTimeSEM._dual_restart!(da, restarted_eps)
    @test ContinuousTimeSEM._dual_final(da) ≈ restarted_eps
    @test ContinuousTimeSEM._dual_final(da) != 1.0
end

@testset "an ordinary (non-restarted) average still converges as before" begin
    # The fix touches only the restart path; an ordinary run of updates should
    # still average as dual averaging does; regression check that nothing
    # about the constructor change broke the common case.
    da = ContinuousTimeSEM._DualAverage(1.0, 0.8)
    local eps
    for _ in 1:200
        eps = ContinuousTimeSEM._dual_update!(da, 0.8)
    end
    # Converged to its own target acceptance with no more forcing, so the
    # final averaged step size should be close to the running one.
    @test ContinuousTimeSEM._dual_final(da) ≈ eps rtol = 0.2
end
