"""
    History{T<:Real,N,M}

History of a `N`-variate point process with marks of type `M` and temporal locations of type `T`.

# Fields

- `times::NTuple{N,Vector{T}}`: each value of the NTuple is a sorted vector of event times
- `marks::NTuple{N,Vector{M}}`: each value of the NTuple is the associated vector of event marks
- `interval::Interval{T}`: observation interval
"""
mutable struct History{T<:Real, N, M}
    times::NTuple{N,Vector{T}}
    marks::NTuple{N,Vector{M}}
    interval::Interval{T}

    function History(times::NTuple{N,Vector{T}}, marks::NTuple{N,Vector{M}}, interval::Interval{T}) where {T<:Real,N,M}
        # Check that times and marks have same length
        if !mapreduce((x,y) -> length(x)==length(y), &, times, marks)
            throw(DimensionMismatch("Times and marks must have same length"))
        end
        # Check that the times are sorted
        if !mapreduce(issorted, &, times)
            @warn "Times have been sorted."
            if N==1
                times, marks = sorttimesandmarks(times[1], marks[1])
                times = (times,)
                marks = (marks,)
            else
                times, marks = map(sorttimesandmarks, times, marks)
            end
        end
        # Check that the times are included in interval
        if !mapreduce(x -> first(x) in interval, &, times)
            throw(ArgumentError("There is at least one time too small: times must be in the interval "*string(interval)*"."))
        end
        if !mapreduce(x -> last(x) in interval, &, times)
            throw(ArgumentError("There is at least one time too large: times must be in the interval "*string(interval)*"."))
        end
        return new{T, N, M}(times, marks, interval)
    end
end

# Constructors
function History(times::NTuple{N,Vector{T}}, marks::NTuple{N,Vector{M}}, interval::Interval{S}) where {T<:Real,S<:Real,N,M}
    NewT = promote_type(T,S)
    times = map(x -> convert(Vector{NewT}, x), times)
    interval = convert(Interval{NewT}, interval)
    return History(times, marks, interval)
end
History(times::Vector{T}, marks::Vector{M}, interval::Interval{S}) where {T<:Real, S<:Real, M} = History((times,), (marks,), interval)
History(times::NTuple{N,Vector{T}}, interval::Interval{S}) where {T<:Real, S<:Real, N} = History(times, ntuple(i -> fill(nothing,length(times[i])), N), interval)
History(times::Vector{T}, interval::Interval{S}) where {T<:Real, S<:Real} = History((times,), interval)

function Base.show(io::IO, h::History{T,N,M}) where {T,N,M}
    N == 1 ? type = "Univariate " : type = string(N)*"-variate "
    M == Nothing ? marked1 = "" : marked1 = "Marked "
    M == Nothing ? marked2 = "" : marked2 = " and marks of type $M"
    print(io, marked1*type*"History with times of type $T"*marked2*" observed on interval ")
    print(io, h.interval)
    print(io, ".")
end

"""
    event_times(h, n)

Return the vector of the event times corresponding to the `n`-th dimension of `h`.
"""
event_times(h::History, n::Integer) = h.times[n]

"""
    event_marks(h, n)

Return the vector of the marks corresponding to the `n`-th dimension of `h` (sorted according to their event times).
"""
event_marks(h::History, n::Integer) = h.marks[n]

"""
    min_time(h)

Return the starting time of `h` (not the same as the first event time).
"""
min_time(h::History) = first(h.interval)

"""
    max_time(h)

Return the end time of `h` (not the same as the last event time).
"""
max_time(h::History) = last(h.interval)

"""
    length(h::History)

Count events in `h`.
"""
Base.length(h::History) = reduce(sum,map(length,h.times))

"""
    length(h, n)

Count events corresponding to the `n`-th dimension of `h`.
"""
Base.size(h::History) = map(length, h.times)

"""
    ndims(h)

Return the dimension of the uni- or multivariate history.
"""
Base.ndims(h::History) = length(h.times)

# TO DO
# Base.eltype()

function intersect_timesandmarks(times::Vector{T},marks::Vector,interval::Interval{T}) where {T<:Real} 
    i_min = searchsortedfirst(times, first(interval))
    i_max = searchsortedlast(times, last(interval))
    return times[i_min:i_max], marks[i_min:i_max]
end

function intersect_times(times::Vector{T},interval::Interval{T}) where {T<:Real} 
    i_min = searchsortedfirst(times, first(interval))
    i_max = searchsortedlast(times, last(interval))
    return times[i_min:i_max]
end

function Base.intersect(h::History{T,N,M},interval::Interval{T}) where {T<:Real,N,M}
    if M == Nothing
        times = map(intersect_times, h.times, ntuple(i -> interval, dimension(h)))
        return History(times, interval)
    else
        timesandmarks = map(intersect_timesandmarks, h.times, h.marks, ntuple(i -> interval, dimension(h)))
        return History(ntuple(i -> timesandmarks[i][1], N), ntuple(i -> timesandmarks[i][2], N), interval)
    end
end


function sorttimesandmarks(times::Vector, marks::Vector)
    order = sortperm(times)
    return times[order], marks[order]
end

# """
#     nb_events(h, tmin, tmax)

# Count events in `h` during the interval `[tmin, tmax)`.
# """
# function nb_events(h::History{M,T}, tmin, tmax) where {M,T}
#     i_min = searchsortedfirst(event_times(h), tmin)
#     i_max = searchsortedlast(event_times(h), tmax - eps(tmax))
#     return i_max - i_min + 1
# end

# """
#     has_events(h)

# Check the presence of events in `h`.
# """
# has_events(h::History) = nb_events(h) > 0

# """
#     has_events(h, tmin, tmax)

# Check the presence of events in `h` during the interval `[tmin, tmax)`.
# """
# has_events(h::History, tmin, tmax) = nb_events(h, tmin, tmax) > 0

# """
#     duration(h)

# Compute the difference `h.tmax - h.tmin`.
# """
# duration(h::History) = max_time(h) - min_time(h)

# """
#     push!(h, t, m)

# Add event `(t, m)` at the end of history `h`.
# """
# function Base.push!(h::History, t, m; check=true)
#     if check
#         @assert h.tmin <= t < h.tmax
#         @assert (length(h) == 0) || (h.times[end] <= t)
#     end
#     push!(h.times, t)
#     push!(h.marks, m)
#     return nothing
# end

# """
#     append!(h1, h2)

# Add all the events of `h2` at the end of `h1`.
# """
# function Base.append!(h1::History, h2::History)
#     max_time(h1) ≈ min_time(h2) || return false
#     append!(h1.times, h2.times)
#     append!(h1.marks, h2.marks)
#     h1.tmax = h2.tmax
#     return true
# end

# """
#     time_change(h, Λ)

# Apply the time rescaling `t -> Λ(t)` to history `h`.
# """
# function time_change(h::History, Λ)
#     new_times = Λ.(event_times(h))
#     new_marks = copy(event_marks(h))
#     new_tmin = Λ(min_time(h))
#     new_tmax = Λ(max_time(h))
#     return History(; times=new_times, marks=new_marks, tmin=new_tmin, tmax=new_tmax)
# end

# """
#     split_into_chunks(h, chunk_duration)

# Split `h` into a vector of consecutive histories with individual duration `chunk_duration`.
# """
# function split_into_chunks(h::History{M}, chunk_duration) where {M}
#     chunks = History{M}[]
#     limits = collect(min_time(h):chunk_duration:max_time(h))
#     if !(limits[end] ≈ max_time(h))
#         push!(limits, max_time(h))
#     end
#     for (a, b) in zip(limits[1:(end - 1)], limits[2:end])
#         times = [t for t in event_times(h) if a <= t < b]
#         marks = [m for (t, m) in zip(event_times(h), event_marks(h)) if a <= t < b]
#         chunk = History(; times=times, marks=marks, tmin=a, tmax=b)
#         push!(chunks, chunk)
#     end
#     return chunks
# end

# """
#     max_mark(h; [init])

# Return the largest event mark if it is larger than `init`, and `init` otherwise.
# """
# max_mark(h::History; init=first(event_marks(h))) = maximum(event_marks(h); init=init)

# """
#     min_mark(h; [init])

# Return the smallest event mark if it is smaller than `init`, and `init` otherwise.
# """
# min_mark(h::History; init=first(event_marks(h))) = minimum(event_marks(h); init=init)
