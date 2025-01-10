function Base.convert(::Type{Interval{T}}, interval::Interval{S,L,R}) where {T, S, L <: Bounded, R <: Bounded}
    Interval{T,L,R}(convert(T, first(interval)), convert(T, last(interval)))
end
function Base.convert(::Type{Interval{T}}, interval::Interval{S,L,R}) where {T, S, L <: Bounded, R <: Unbounded}
    Interval{T,L,R}(convert(T, first(interval)), nothing)
end
function Base.convert(::Type{Interval{T}}, interval::Interval{S,L,R}) where {T, S, L <: Unbounded, R <: Bounded}
    Interval{T,L,R}(nothing, convert(T, last(interval)))
end
