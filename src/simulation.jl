"""
    simulate_ogata(rng, pp, tmin, tmax)

Simulate a temporal point process `pp` on interval `[tmin, tmax)` using Ogata's algorithm.
"""
function simulate_ogata(
    rng::AbstractRNG, pp::AbstractPointProcess, tmin::T, tmax::T
) where {T<:Real}
    h = History(; times=T[], marks=[], tmin=tmin, tmax=tmax)
    M = typeof(rand(mark_distribution(pp, zero(T), h)))
    h.marks = M[]
    t = tmin
    while t < tmax
        B, L = ground_intensity_bound(pp, t + eps(t), h)
        τ = B > 0 ? rand(rng, Exponential(inv(B))) : typemax(inv(B))
        if τ > L
            t = t + L
        elseif τ <= L
            U_max = ground_intensity(pp, t + τ, h) / B
            U = rand(rng, typeof(U_max))
            if U < U_max
                m = rand(rng, mark_distribution(pp, t + τ, h))
                if t + τ < tmax
                    push!(h, t + τ, m)
                end
            end
            t = t + τ
        end
    end
    return h
end

"""
    rand([rng,] pp, tmin, tmax)

Alias for `simulate_ogata`.
"""
function Base.rand(rng::AbstractRNG, pp::AbstractPointProcess, tmin, tmax)
    return simulate_ogata(rng, pp, tmin, tmax)
end

function Base.rand(pp::AbstractPointProcess, args...; kwargs...)
    return rand(default_rng(), pp, args...; kwargs...)
end
