"""
    MarkovModulatedPoissonProcess{M,Tr<:ContinuousMarkovChain,Em<:PoissonProcess{M}}

Markov-Modulated Poisson Process with mark type `M`.

# Fields
- `transitions::Tr`: state evolution process.
- `emissions::Vector{Em}`: one emission distribution per state.
"""
@with_kw struct MarkovModulatedPoissonProcess{
    M,
    Tr<:ContinuousMarkovChain,
    Em<:PoissonProcess{M},
}
    transitions::Tr
    emissions::Vector{Em}
end

## Access

transitions(mmpp::MarkovModulatedPoissonProcess) = mmpp.transitions
initial_distribution(mmpp::MarkovModulatedPoissonProcess) =
    initial_distribution(transitions(mmpp))
transition_matrix(mmpp::MarkovModulatedPoissonProcess) =
    transition_matrix(transitions(mmpp))
rate_matrix(mmpp::MarkovModulatedPoissonProcess) = rate_matrix(transitions(mmpp))

emissions(mmpp::MarkovModulatedPoissonProcess) = mmpp.emissions
emission(mmpp::MarkovModulatedPoissonProcess, s::Int) = mmpp.emissions[s]
nstates(mmpp::MarkovModulatedPoissonProcess) = length(emissions(mmpp))

## Simulation

function Base.rand(
    rng::AbstractRNG,
    mmpp::MarkovModulatedPoissonProcess{M},
    tmin,
    tmax,
) where {M}
    state_history = rand(transitions(mmpp), tmin, tmax)
    transition_times, states = event_times(state_history), event_marks(state_history)
    observations = TemporalHistory(times = Float64[], marks = M[], tmin = tmin, tmax = tmax)
    for k = 1:length(transition_times)
        local_s = states[k]
        local_tmin = transition_times[k]
        local_tmax = k < length(transition_times) ? transition_times[k+1] : tmax
        local_observations = rand(emission(mmpp, local_s), local_tmin, local_tmax)
        append!(observations, local_observations)
    end
    return state_history, observations
end

Base.rand(mmpp::MarkovModulatedPoissonProcess, args...) =
    rand(Random.GLOBAL_RNG, mmpp, args...)
