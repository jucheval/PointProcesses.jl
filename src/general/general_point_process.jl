struct PointProcess{M} <: AbstractPoissonProcess{M}
    ground_intensity::Function
    ground_intensity_bound::Function
    mark_distribution::Function
end

# Constructor
function PointProcess(ground_intensity, ground_intensity_bound, mark_distribution)
    arbitrary_mark_distribution = mark_distribution(0.0, History(Float64[], [], 0.0, 1.0))
    M = eltype(arbitrary_mark_distribution)
    PointProcess{M}(ground_intensity, ground_intensity_bound, mark_distribution)
end

# Show
function Base.show(io::IO, pp::PointProcess{M})
    return print(io, "Generic point process with marks of type $M")
end

# Intensity functions
ground_intensity(pp::PointProcess, t, h) = pp.ground_intensity(t, h)
mark_distribution(pp::PointProcess, t, h) = pp.mark_distribution(t, h)
ground_intensity_bound(pp::PointProcess, t, h) = pp.ground_intensity_bound(t, h)
integrated_ground_intensity(pp::PointProcess, h, a, b) = error("No generic method to compute the exact integrated ground intensity")

# Convert functions
Base.convert(PointProcess, pp::AbstractPointProcess) = PointProcess(ground_intensity)