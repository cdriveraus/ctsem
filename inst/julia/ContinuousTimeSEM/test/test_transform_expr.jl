# Transforms held as data (src/transform_expr.jl): an expression must give
# exactly what the closure it replaces gives, value and seeded derivative, and
# models that differ only in their expressions must be one type.
const _TX = ContinuousTimeSEM
const _TXD = ForwardDiff.Dual{Nothing,Float64,1}
_tx_dual(x, p) = _TXD(x, ForwardDiff.Partials((p,)))
# A closure `eval`ed while this file runs is newer than the running code, so it
# is called through `invokelatest`.
_tx_call(f, c) = Base.invokelatest(f, c)
# Value and seeded partial, exactly. A closure reading nothing seeded returns a
# plain real where an expression returns the context's dual type with a zero
# partial; the reverse pass reads both the same way (`_partial1_as`).
_tx_same(a, b) = ForwardDiff.value(a) === ForwardDiff.value(b) &&
    _TX._partial1(a) === _TX._partial1(b)

@testset "CTSEMExpr evaluates as its closure does" begin
    matrix = ["DRIFT", "DRIFT", "DRIFT", "DRIFT", "PARS", "PARS", "LAMBDA", "LAMBDA"]
    row = [1, 2, 1, 2, 1, 2, 1, 1]
    col = [1, 1, 2, 2, 1, 1, 1, 2]
    layout = _TX.CTSEMLayout(_TX.retrieve_layout(matrix, row, col)...)
    data = collect(range(0.3, step = 0.11, length = layout.total))
    state = [0.7, -0.4, 1.3]
    tdpreds = [0.25, -1.5]
    tipreds = [2.0]
    ctx = _TX.CTSEMRowContext(state, _TX.CTSEMPars(data, layout), tdpreds, tipreds,
        0.5, 0.25, 1, 3)
    sources = [
        "10*state[3]",
        "state[1] + PARS[1]",
        "PARS[1,1]*(1+0.5*state[1])",
        "DRIFT[2,1]*state[2]^2 - exp(PARS[2,1])/state[1]",
        "-state[2]*ctx.tdpreds[2]",
        "log1p_exp(state[1]) + inv_logit(PARS[2]) * ctx.tipreds[1]",
        "sqrt(abs(state[2])) + tanh(LAMBDA[1,2]) + state[1]^-1 + 2^3*state[3]",
        "max(state[1], PARS[1]) + min(state[2], 0.1) + state[1] * state[2] * state[3]",
        "(state[1] + state[2] + state[3]) / 3",
        "state[1]^0.5 * log(DRIFT[1,2]) - expm1(state[2]) + cos(state[3])^2",
        "state[1]/3 + PARS[1]*ctx.tdpreds[1]/4 - (ctx.tdpreds[2]*2 + 1)*state[2]",
        "exp(ctx.tipreds[1]/2)*PARS[2,1] + 1/(1 + exp(-state[3]))",
        # Past eight nodes: evaluated in the task's buffers.
        "-PARS[1,1]*(1 + 0.5*state[1]^2 - 0.2*state[1]^3 + 0.1*state[1]^4)/(1 + exp(-PARS[2,1]*state[1])) - LAMBDA[1,1]*sin(state[2])^2*ctx.tdpreds[1]",
    ]
    for s in sources
        e = _TX._ctsem_expr(s, layout)
        @test e isa _TX.CTSEMExpr
        f = _TX._transform_closure(_TX.generate_complex_transform_string(s))
        @test e(ctx) === _tx_call(f, ctx)
        # Seeded on each state and on each parameter cell in turn, as the
        # reverse pass differentiates a transform: state and parameters in
        # one dual type, as `_sync_dual_context!` builds them.
        for k in eachindex(state)
            dstate = _tx_dual.(state, 0.0)
            dstate[k] = _tx_dual(state[k], 1.0)
            dctx = _TX.CTSEMRowContext(dstate, _TX.CTSEMPars(_tx_dual.(data, 0.0), layout),
                tdpreds, tipreds, 0.5, 0.25, 1, 3)
            @test _tx_same(e(dctx), _tx_call(f, dctx))
        end
        for j in eachindex(data)
            ddata = _tx_dual.(data, 0.0)
            ddata[j] = _tx_dual(data[j], 1.0)
            dctx = _TX.CTSEMRowContext(_tx_dual.(state, 0.0), _TX.CTSEMPars(ddata, layout),
                tdpreds, tipreds, 0.5, 0.25, 1, 3)
            @test _tx_same(e(dctx), _tx_call(f, dctx))
        end
    end
    # Outside the vocabulary, or a cell the model does not have: the closure.
    @test _TX._ctsem_expr("ifelse(state[1] > 0, 1.0, 2.0)", layout) === nothing
    @test _TX._ctsem_expr("TDPREDEFFECT[1,1]*state[1]", layout) === nothing
    @test _TX._ctsem_expr("DRIFT[3,1]*state[1]", layout) === nothing
    @test _TX._ctsem_expr("state[i]", layout) === nothing
    # Without a layout, only whether it is in the vocabulary.
    @test _TX._ctsem_expr("10*state[3] + TDPREDEFFECT[1,1]", nothing) isa _TX.CTSEMExpr
end

@testset "An expression's pullback gives the seeded duals' partials" begin
    # The reverse pass takes every partial of an expression in one sweep
    # (`_ctsem_expr_pullback!`) where the closures took one seeded evaluation
    # per input; the two agree to rounding. A random effect's `c*state[k]` is
    # not swept but read off, and is the dual's to the bit.
    matrix = ["DRIFT", "DRIFT", "DRIFT", "DRIFT", "PARS", "PARS", "LAMBDA", "LAMBDA"]
    row = [1, 2, 1, 2, 1, 2, 1, 1]
    col = [1, 1, 2, 2, 1, 1, 1, 2]
    layout = _TX.CTSEMLayout(_TX.retrieve_layout(matrix, row, col)...)
    data = collect(range(0.3, step = 0.11, length = layout.total))
    state = [0.7, -0.4, 1.3]
    tdpreds = [0.25, -1.5]
    tipreds = [2.0]
    zd(v) = _tx_dual.(v, 0.0)
    dctx = _TX.CTSEMRowContext(zd(state), _TX.CTSEMPars(zd(data), layout), tdpreds,
        zd(tipreds), 0.5, 0.25, 1, 3)
    sources = [
        "PARS[1,1]*(1+0.5*state[1])",
        "DRIFT[2,1]*state[2]^2 - exp(PARS[2,1])/state[1]",
        "-state[2]*ctx.tdpreds[2]",
        "log1p_exp(state[1]) + inv_logit(PARS[2]) * ctx.tipreds[1]",
        "sqrt(abs(state[2])) + tanh(LAMBDA[1,2]) + state[1]^-1",
        "max(state[1], PARS[1]) + min(state[2], 0.1)",
        "state[1]^0.5 * log(DRIFT[1,2]) - expm1(state[2])",
        "exp(ctx.tipreds[1]/2)*PARS[2,1] + 1/(1 + exp(-state[3]))",
        "cos(state[3])^2 - sin(state[1]) * logit(0.3 + 0.1*state[1])",
        "-PARS[1,1]*(1 + 0.5*state[1]^2 - 0.2*state[1]^3 + 0.1*state[1]^4)/(1 + exp(-PARS[2,1]*state[1])) - LAMBDA[1,1]*sin(state[2])^2*ctx.tdpreds[1]",
    ]
    seeded(e, k, which) = begin
        d = which === :state ? (s = zd(state); s[k] = _tx_dual(state[k], 1.0);
                _TX.CTSEMRowContext(s, _TX.CTSEMPars(zd(data), layout), tdpreds, zd(tipreds), 0.5, 0.25, 1, 3)) :
            which === :par ? (p = zd(data); p[k] = _tx_dual(data[k], 1.0);
                _TX.CTSEMRowContext(zd(state), _TX.CTSEMPars(p, layout), tdpreds, zd(tipreds), 0.5, 0.25, 1, 3)) :
            (t = zd(tipreds); t[k] = _tx_dual(tipreds[k], 1.0);
                _TX.CTSEMRowContext(zd(state), _TX.CTSEMPars(zd(data), layout), tdpreds, t, 0.5, 0.25, 1, 3))
        _TX._partial1(e(d))
    end
    for s in sources
        e = _TX._ctsem_expr(s, layout)
        pbar = zeros(length(data))
        sbar = zeros(length(state))
        ti = (supports = Int[1], bar = zeros(1), sampled = [true])
        n = length(e.op)
        g = zeros(n)
        g[e.root] = 1.0
        _TX._ctsem_expr_pullback!(pbar, sbar, e, dctx, ti, 1.0, Float64, zeros(n), g)
        for j in eachindex(data)
            @test pbar[j] ≈ seeded(e, j, :par) rtol = 1e-13 atol = 1e-15
        end
        for k in eachindex(state)
            @test sbar[k] ≈ seeded(e, k, :state) rtol = 1e-13 atol = 1e-15
        end
        @test ti.bar[1] ≈ seeded(e, 1, :ti) rtol = 1e-13 atol = 1e-15
    end
    carrier = _TX._ctsem_expr("10*state[3]", layout)
    @test carrier.form === _TX._XFORM_SCALED_STATE
    @test carrier.fc * 1.0 === seeded(carrier, 3, :state)
end

@testset "Regular transforms outside the default templates are expressions" begin
    values = [0.3, -0.7, 1.1]
    texts = [
        "10 * log1p_exp(2 * param[3]) + 1e-10",    # a default template
        "-(1e-06 + 2 * log1p_exp(-(2 * param[2])))", # a default template
        "3 * log1p_exp(2 * param[2]) + 1e-10",     # a scale that is not
        "exp(param[1]) * 0.5",
        "1/(1 + exp(-(4 * param[3] - 1)))",
        "param[2] + 0.25 * param[3]",              # composed from two
    ]
    for t in texts
        g = _TX._regular_transform_closure(t)
        default = _TX._regular_transform_index(t) !== nothing && _TX._ctsem_default_template(t)
        @test (g isa _TX.CTSEMExpr) == !default
        f = _TX._transform_closure(_TX.generate_transform_string(t))
        @test _tx_call(g, values) === _tx_call(f, values)
        for k in eachindex(values)
            dv = _tx_dual.(values, 0.0)
            dv[k] = _tx_dual(values[k], 1.0)
            @test _tx_same(_tx_call(g, dv), _tx_call(f, dv))
        end
    end
    # A transform reading no parameter keeps its closure, so the support probe
    # reports it as the rendering fault it is.
    @test _TX._ctsem_regular_expr("2.5") === nothing
end

@testset "Models differing only in their expressions are one type" begin
    base = (matrix = [:A, :A, :PARS], row = [1, 2, 1], col = [1, 1, 1],
        parnumber = Union{Missing,Int}[1, missing, 2],
        value = Union{Missing,Float64}[missing, missing, missing])
    default_tf = Union{Missing,String}["param[1]", missing, "param[2]"]
    one = DataFrame(; base..., transform = default_tf,
        predicttransform = Union{Missing,String}[missing, "state[1] + PARS[1]", missing],
        updatetransform = Union{Missing,String}[missing, missing, missing])
    two = DataFrame(; base..., transform = default_tf,
        predicttransform = Union{Missing,String}[missing, "10*state[2]", missing],
        updatetransform = Union{Missing,String}[missing, "PARS[1]*state[1]^2", missing])
    none = DataFrame(; base..., transform = default_tf,
        predicttransform = Union{Missing,String}[missing, missing, missing],
        updatetransform = Union{Missing,String}[missing, missing, missing])
    # Parameter transforms of the user's own, outside the default templates.
    custom = DataFrame(; base...,
        transform = Union{Missing,String}["3 * log1p_exp(param[1])", missing, "exp(param[2]) - 1"],
        predicttransform = Union{Missing,String}[missing, missing, missing],
        updatetransform = Union{Missing,String}[missing, missing, missing])
    specs = [ekf_from_data_frame(df) for df in (one, two, none, custom)]
    @test specs[1].predict_transforms isa Vector{_TX.CTSEMExpr}
    @test isempty(specs[3].predict_transforms)
    @test allequal(typeof.(specs))
    # Compiled on request: closures, a type of the model's own, the same
    # values, and cells read by flat index (`_ctsem_flat_cells`).
    compiled = ContinuousTimeSEM.ekf_from_columns(two.matrix, two.row, two.col,
        _tbl_int(two.parnumber), _tbl_float(two.value), _tbl_string(two.transform),
        _tbl_string(two.predicttransform), _tbl_string(two.updatetransform),
        fill("", nrow(two)); compile_transforms = true)
    @test compiled.predict_transforms isa Tuple
    @test typeof(compiled) != typeof(specs[2])
    layout = specs[2].parameter_axis
    @test _TX._ctsem_flat_cells("ctx.pars.PARS[1,1] * ctx.state[1]", layout) ==
        "ctx.pars[$(layout.first[_TX._layout_index(layout, :PARS)])] * ctx.state[1]"
    @test _TX._ctsem_flat_cells("ctx.pars.NOSUCH[1,1]", layout) == "ctx.pars.NOSUCH[1,1]"
    vals = collect(range(0.4, step = 0.1, length = layout.total))
    ctx = _TX.CTSEMRowContext([0.7, -0.3], _TX.CTSEMPars(vals, layout), Float64[], Float64[],
        0.5, 0.25, 1, 2)
    for (e, f) in zip(Iterators.flatten((specs[2].predict_transforms, specs[2].update_transforms)),
                      Iterators.flatten((compiled.predict_transforms, compiled.update_transforms)))
        @test e(ctx) === _tx_call(f, ctx)
    end
    # A transform outside the vocabulary keeps its closure, and its own type.
    odd = DataFrame(; base..., transform = default_tf,
        predicttransform = Union{Missing,String}[missing, "ifelse(state[1] > 0, 1.0, 2.0)", missing],
        updatetransform = Union{Missing,String}[missing, missing, missing])
    @test ekf_from_data_frame(odd).predict_transforms isa Tuple
end
