## A content-keyed cache in front of the covariance construction, for BOTH
## routes.
##
## `_ekf_predict_step!` and its neighbours call `sdcovsqrt2cov!` unconditionally
## on every row, and the julia filter has no equivalent of the stan writer's
## `whenmat`/`statedep` gate -- so T0VAR, DIFFUSION and MANIFESTVAR are rebuilt
## whether or not anything about them changed. Measured on a k=3, 30-subject,
## 300-row model with individual differences in one DIFFUSION correlation: 6000
## construction calls per gradient on the laplace route for 140 distinct
## matrices.
##
## Keyed on `mat`'s full contents rather than a hash, so an entry cannot go
## stale. A hit costs up to `_COVCACHE_SLOTS` comparisons of O(d^2); the
## construction it replaces is O(d^3). A table rather than one slot because the
## three matrices of a given size share the entry and would otherwise evict each
## other.
##
## Disabled by `ctsem_cov_cache!(false)`, which exists so a test can measure the
## uncached cost rather than infer it.

const _CTSEM_COV_CACHE = Ref(true)

"""Enable or disable the covariance construction cache. Returns the previous setting."""
function ctsem_cov_cache!(on::Bool)
    prev = _CTSEM_COV_CACHE[]
    _CTSEM_COV_CACHE[] = on
    return prev
end
# Bang-free aliases, as JuliaConnectoR addresses module members by name.
ctsem_cov_cache(on::Bool) = ctsem_cov_cache!(on)
ctsem_cov_cache() = _CTSEM_COV_CACHE[]

const _COVCACHE_SLOTS = 8
const _COVCACHE_CALLS = Ref(0)
const _COVCACHE_MISSES = Ref(0)

"""Calls into the covariance construction, and how many missed the cache."""
ctsem_cov_cache_counts() = (calls = _COVCACHE_CALLS[], misses = _COVCACHE_MISSES[])

"""Zero the covariance construction counters."""
function ctsem_cov_cache_reset_counts()
    _COVCACHE_CALLS[] = 0
    _COVCACHE_MISSES[] = 0
    return nothing
end

struct CovCache{T,D}
    keys::Vector{Matrix{T}}
    outs::Vector{Matrix{T}}
    nslots::Base.RefValue{Int}
    nextslot::Base.RefValue{Int}
end

CovCache{T,D}() where {T,D} = CovCache{T,D}(
    [zeros(T, D, D) for _ in 1:_COVCACHE_SLOTS],
    [zeros(T, D, D) for _ in 1:_COVCACHE_SLOTS], Ref(0), Ref(1))

const _COVCACHE = Dict{Tuple{DataType,Int,Int,Bool},Any}()
const _COVCACHE_LOCK = ReentrantLock()

# The construction in force is part of the key. Without it an entry built under
# one route is served under the other, which a switch mid-session -- every
# comparison test and benchmark here -- turns into a wrong answer with no
# symptom. Caught by the transparency test, not by inspection.
function _covcache(::Type{T}, ::Val{d}) where {T,d}
    key = (T, d, Threads.threadid(), _CTSEM_COV_EXPM[])
    c = get(_COVCACHE, key, nothing)
    if c === nothing
        lock(_COVCACHE_LOCK) do
            c = get(_COVCACHE, key, nothing)
            if c === nothing
                c = CovCache{T,d}()
                _COVCACHE[key] = c
            end
        end
    end
    return c::CovCache{T,d}
end

# The whole of `mat` is the key: the diagonal carries the standard deviations
# and the lower triangle the off-diagonal coordinates, and both feed the result.
@inline function _covcache_lookup(c::CovCache{T,D}, mat, ::Val{d}) where {T,D,d}
    @inbounds for e in 1:c.nslots[]
        k = c.keys[e]
        same = true
        for j in 1:d, i in 1:d
            if k[i, j] != mat[i, j]
                same = false
                break
            end
        end
        same && return e
    end
    return 0
end

@inline function _covcache_store!(c::CovCache{T,D}, mat, out, ::Val{d}) where {T,D,d}
    slot = c.nextslot[]
    @inbounds begin
        for j in 1:d, i in 1:d
            c.keys[slot][i, j] = mat[i, j]
        end
        copyto!(c.outs[slot], out)
    end
    c.nslots[] = max(c.nslots[], slot)
    c.nextslot[] = slot == _COVCACHE_SLOTS ? 1 : slot + 1
    return nothing
end
