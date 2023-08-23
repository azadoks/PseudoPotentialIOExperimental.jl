abstract type RadialMesh{T} <: AbstractVector{T} end

function RadialMesh(x::AbstractVector{T}) where {T}
    for mesh_type in (UniformMesh, LogMesh)
        mesh = mesh_type(x)
        all(mesh .≈ x) && return mesh
    end

    error("Could not construct a known mesh type which matches input data.")
end

function RadialMesh(x::AbstractVector{T}, dx::AbstractVector{T}) where {T}
    for mesh_type in (UniformMesh, LogMesh)
        mesh = mesh_type(x)
        all(mesh .≈ x) && all(deriv(mesh) .≈ dx) && return mesh
    end

    error("Could not construct a known mesh type which matches input data.")
end

function RadialMesh(x::StepRangeLen)
    return UniformMesh(x)
end

Base.length(mesh::RadialMesh) = mesh.n
Base.size(mesh::RadialMesh) = (mesh.n,)
Base.firstindex(mesh::RadialMesh) = 1
Base.lastindex(mesh::RadialMesh) = mesh.n
Base.similar(mesh::RadialMesh{T}) where {T} = Vector{T}(undef, mesh.n)
Base.similar(mesh::RadialMesh, ::Type{S}) where {S} = Vector{S}(undef, mesh.n)
Base.eachindex(mesh::RadialMesh) = Base.OneTo(mesh.n)
Base.iterate(mesh::RadialMesh, state=1) = state > mesh.n ? nothing : (mesh[state], state+1)
Base.eltype(::Type{<:RadialMesh{T}}) where {T} = T
Base.summary(io::IO, mesh::RadialMesh{T}) where {T} = print(io, "$(mesh.n)-element $(typeof(mesh))")

function Base.getindex(mesh::RadialMesh{T}, i::AbstractRange)::Vector{T} where {T}
    map(Base.Fix1(getindex, mesh), i)
end

function Base.getindex(mesh::RadialMesh{T}, ::Colon)::Vector{T} where{T}
    Base.getindex(mesh, eachindex(mesh))
end

function (mesh::RadialMesh{T})(i::Integer)::T where {T}
    return mesh[i]
end

function Base.diff(mesh::RadialMesh{T}) where {T}
    return map(Base.Fix1(diff, mesh), 1:mesh.n-1)
end

function Base.diff(mesh::RadialMesh{T}, i::AbstractUnitRange) where {T}
    return map(Base.Fix1(diff, mesh), i)
end

function deriv(mesh::RadialMesh{T}) where {T}
    return map(Base.Fix1(deriv, mesh), eachindex(mesh))
end

function deriv(mesh::RadialMesh{T}, i::AbstractUnitRange) where {T}
    return map(Base.Fix1(deriv, mesh), i)
end

struct ArbitraryMesh{T} <: RadialMesh{T}
    x::AbstractVector{T}
    dx::AbstractVector{T}
    Δx::AbstractVector{T}
    n::Int
end

function ArbitraryMesh(x::AbstractVector, dx::AbstractVector)
    length(x) == length(dx) || throw(ArgumentError("x and dx must contain the same number of elements"))

    Δx = diff(x)
    n = length(x)
    return ArbitraryMesh(x, dx, Δx, n)
end

function ArbitraryMesh(x::AbstractVector)
    fit = RadialMesh(x)
    dx = deriv(fit)
    Δx = diff(x)
    n = length(x)
    return ArbitraryMesh(x, dx, Δx, n)
end

function Base.getindex(mesh::ArbitraryMesh{T}, i::Integer)::T where {T}
    return mesh.x[i]
end

function Base.diff(mesh::ArbitraryMesh{T}, i::Integer)::T where {T}
    return mesh.Δx[i]
end

function deriv(mesh::ArbitraryMesh{T}, i::Integer)::T where {T}
    return mesh.dx[i]
end

function index_leq(mesh::ArbitraryMesh{T}, x::T)::Int where {T}
    return findfirst(≤(x), mesh.x)
end

# # r(i) = b + a * (i - 1)
# # Δr(i) = r(i + 1) - r(i) = [b + a * (i + 1 - 1)] - [b + a * (i - 1)] = a
# # dr(i)/di = i
struct UniformMesh{T} <: RadialMesh{T}
    a::T
    b::T
    x1::T
    xn::T
    n::Int
end

Base.show(io::IO, mesh::UniformMesh) = print(io, "$(mesh.x1):$(mesh.a):$(mesh.xn)")

function UniformMesh(x::AbstractVector)
    return UniformMesh(first(x), last(x), length(x))
end

function UniformMesh(x1::T, xn::T, maximal_spacing::T) where {T<:Real}
    n = ceil(Int, (xn - x1) / maximal_spacing)
    return UniformMesh(x1, xn, n)
end

function UniformMesh(x1::T, xn::T, n::Integer) where {T<:Real}
    n > 1 || throw(ArgumentError("n must be greater than 1"))
    xn > x1 || throw(ArgumentError("xn must be strictly greater than x1"))
    a = (xn - x1) / (n - 1)
    b = x1
    M = promote_type(typeof(a), typeof(b))
    return UniformMesh{M}(a, b, x1, xn, n)
end

function Base.getindex(mesh::UniformMesh{T}, i::Integer)::T where {T}
    # This check  is important but _very_ slow
    # 1 < i <= mesh.n + 1 || throw(BoundsError(mesh, i))
    return mesh.b + mesh.a * (i - 1)
end

function Base.diff(mesh::UniformMesh{T}, i::Integer)::T where {T}
    0 < i < mesh.n || throw(BoundsError(mesh, i+1))
    return mesh.a
end

function deriv(mesh::UniformMesh{T}, i::Integer)::T where {T}
    1 <= i <= mesh.n || throw(BoundsError(mesh, i))
    return mesh.a
end

function index_leq(mesh::UniformMesh{T}, x::T)::Int where {T}
    return floor(Int, (x - mesh.b) / mesh.a)
end

abstract type LogMesh{T} <: RadialMesh{T} end

function LogMesh(x::AbstractVector)
    iszero(first(x)) ? LogMeshWithZero(x) : LogMeshWithoutZero(x)
end


# r(i) = b exp(a * (i - 1))
# Δr(i) = r(i + 1) - r(i) = b exp(a * i) - b exp(a * (i - 1))
#       = b (exp(a * i) - exp(a * (i - 1)))
# dr(i)/di = a * b * exp(a * (i - 1))
struct LogMeshWithoutZero{T} <: LogMesh{T}
    a::T
    b::T
    x1::T
    xn::T
    n::Int
end

Base.show(io::IO, mesh::LogMeshWithoutZero) = print(io, "x(i) = $(mesh.b) exp( $(mesh.a) (i - 1) )")

function LogMeshWithoutZero(x::AbstractVector)
    return LogMeshWithoutZero(first(x), last(x), length(x))
end

function LogMeshWithoutZero(x1::T, xn::T, n::Integer) where {T}
    n >= 2 || throw(ArgumentError("n must be at least 2"))
    xn > x1 || throw(ArgumentError("xn must be strictly greater than x1"))
    x1 > 0.0 || throw(ArgumentError("x1 must be strictly greater than 0"))
    b = x1
    a = log(xn / x1) / (n - 1)
    return LogMeshWithoutZero{T}(a, b, x1, xn, n)
end

function Base.getindex(mesh::LogMeshWithoutZero{T}, i::Integer)::T where {T}
    # This check is important but _very_ slow
    # 1 <= i <= mesh.n || throw(BoundsError(mesh, i))
    return mesh.b * exp(mesh.a * (i - 1))
end

function Base.diff(mesh::LogMeshWithoutZero{T}, i::Integer)::T where {T}
    0 < i < mesh.n || throw(BoundsError(mesh, i+1))
    return mesh.b * (exp(mesh.a * i) - exp(mesh.a * (i - 1)))
end

function deriv(mesh::LogMeshWithoutZero{T}, i::Integer)::T where {T}
    1 <= i <= mesh.n || throw(BoundsError(mesh, i))
    return mesh.a * mesh.b * exp(mesh.a * (i - 1))
end

function index_leq(mesh::LogMeshWithoutZero{T}, x::T)::Integer where {T}
    return floor(Int, log(x / mesh.b) / mesh.a) + 1
end


# r(i) = b [exp(a * (i - 1)) - 1]
# Δr(i) = r(i + 1) - r(i) = b [exp(a * i) - 1] - b [exp(a * (i - 1)) - 1]
#       = b (exp(a * i) - exp(a * (i - 1)))
# dr(i)/di = a * b * exp(a * (i - 1)) + a * b
struct LogMeshWithZero{T} <: LogMesh{T}
    a::T
    b::T
    x2::T
    x3::T
    n::Int
end

function LogMeshWithZero(x::AbstractVector)
    iszero(first(x)) || throw(ArgumentError("The first element of x must be zero"))
    length(x) >= 3 || throw(ArgumentError("x must have at least 3 elements"))
    return LogMeshWithZero(x[begin+1], x[begin+2], length(x))
end

function LogMeshWithZero(x2::T, x3::T, n::Integer) where {T}
    n >= 2 || throw(ArgumentError("n must be at least 2"))
    b = x2^2 / (x3 - 2x2)
    a = log((x2 + b) / b)
    return LogMeshWithZero{T}(a, b, x2, x3, n)
end

Base.show(io::IO, mesh::LogMeshWithZero) = print(io, "x(i) = $(mesh.b) exp( $(mesh.a) (i - 1) ) - 1")

function Base.getindex(mesh::LogMeshWithZero{T}, i::Integer)::T where {T}
    # This check is important but _very_ slow
    # 1 <= i <= mesh.n || throw(BoundsError(mesh, i))
    return mesh.b * (exp(mesh.a * (i - 1)) - 1)
end

function Base.diff(mesh::LogMeshWithZero{T}, i::Integer)::T where {T}
    0 < i < mesh.n || throw(BoundsError(mesh, i+1))
    return mesh.b * (exp(mesh.a * i) - exp(mesh.a * (i - 1)))
end

function deriv(mesh::LogMeshWithZero{T}, i::Integer)::T where {T}
    1 <= i <= mesh.n || throw(BoundsError(mesh, i))
    return mesh.a * mesh.b * (exp(mesh.a * (i - 1)) - 1) + mesh.a * mesh.b
end

function index_leq(mesh::LogMeshWithZero{T}, x::T)::Int where {T}
    return floor(Int, log(x / mesh.b + 1) / mesh.a) + 1
end
