using Infiltrator
import ForwardDiff: Dual

# TODO Based on this type do some generic things like use spin-scaling relations
# Note: Kind is needed because GGA enhancement or spin scaling etc. work differently
#       for exchange and correlation ... if it does not help, remove it again
abstract type Functional{Family,Kind} end

"""Return the family of a functional. Results are `:lda`, `:gga`, `:mgga` and
`:mggal` (Meta-GGA requiring Laplacian of ρ)"""
family(::Functional{F}) where {F} = F

"""
Return the functional kind: `:x` (exchange), `:c` (correlation), `:k` (kinetic) or
`:xc` (exchange and correlation combined)
"""
kind(::Functional{F,K}) where {F,K} = K

"""Return the identifier corresponding to a functional"""
function identifier end
Base.show(io::IO, fun::Functional) = print(io, identifier(fun))

@doc raw"""
True if the functional needs ``σ = 𝛁ρ ⋅ 𝛁ρ``.
"""
needs_σ(::Functional{F}) where {F} = (F in (:gga, :mgga, :mggal))

@doc raw"""
True if the functional needs ``τ`` (kinetic energy density).
"""
needs_τ(::Functional{F}) where {F} = (F in (:mgga, :mggal))

@doc raw"""
True if the functional needs ``Δ ρ``.
"""
needs_Δρ(::Functional{F}) where {F} = (F in (:mggal,))

"""
Does this functional support energy evaluations? Some don't, in which case
energy terms will not be returned by `potential_terms` and `kernel_terms`,
i.e. `e` will be `false` (a strong zero).
"""
has_energy(::Functional) = true

"""
Return adjustable parameters of the functional and their values.
"""
parameters(::Functional) = ComponentArray{Bool}()

"""Return the type used to represent the [`parameters`](@ref) in the functional"""
parameter_type(f::Functional) = eltype(parameters(f))

"""
Return a new version of the passed functional with its parameters adjusted.
This may not be a copy in case no changes are done to its internal parameters.
Generally the identifier of the functional will be changed to reflect the
change in parameter values unless `keep_identifier` is true.
To get the tuple of adjustable parameters and their current values check out
[`parameters`](@ref). It is not checked that the correct parameters are passed.
"""
change_parameters(f::Functional, ::AbstractArray; keep_identifier=false) = f

# TODO These values are read-only for now and their defaults hard-coded for Float64
"""
Threshold for the density (below this value, functionals and derivatives
evaluate to zero). The threshold may depend on the floating-point type used
to represent densities and potentials, which is passed as the second argument.
"""
threshold_ρ(::Functional, T=Float64) = T(1e-15)  # TODO This might differ between functionals
threshold_σ(f::Functional, T=Float64) = threshold_ρ(f, T)^(4 // 3)
threshold_τ(::Functional, T=Float64)  = T(1e-20)
threshold_ζ(::Functional, T=Float64)  = eps(T)

# Drop dual types from threshold functions
threshold_ρ(f::Functional, T::Type{<:Dual}) = threshold_ρ(f, ForwardDiff.valtype(T))
threshold_σ(f::Functional, T::Type{<:Dual}) = threshold_σ(f, ForwardDiff.valtype(T))
threshold_τ(f::Functional, T::Type{<:Dual}) = threshold_τ(f, ForwardDiff.valtype(T))
threshold_ζ(f::Functional, T::Type{<:Dual}) = threshold_ζ(f, ForwardDiff.valtype(T))

# Silently drop extra arguments from evaluation functions
for fun in (:potential_terms, :kernel_terms)
    @eval begin
        $fun(func::Functional{:lda}, ρ, σ, args...)         = $fun(func, ρ)
        $fun(func::Functional{:gga}, ρ, σ, τ, args...)      = $fun(func, ρ, σ)
        $fun(func::Functional{:mgga}, ρ, σ, τ, Δρ, args...) = $fun(func, ρ, σ, τ)
    end
end

@doc raw"""
    potential_terms(f::Functional, ρ, [σ, τ, Δρ])

Evaluate energy and potential terms at a real-space grid of densities, density
derivatives etc. Not required derivatives for the functional type will be ignored.
Returns a named tuple with keys `e` (Energy per unit volume),
`Vρ` (``\frac{∂e}{∂ρ}``), `Vσ` (``\frac{∂e}{∂σ}``),
`Vτ` (``\frac{∂e}{∂τ}``), `Vl` (``\frac{∂e}{∂(Δρ)}``).
"""
function potential_terms end

@doc raw"""
    kernel_terms(f::Functional, ρ, [σ, τ, Δρ])

Evaluate energy, potential and kernel terms at a real-space grid of densities, density
derivatives etc. Not required derivatives for the functional type will be ignored.
Returns a named tuple with the same keys as `potential_terms` and additionally
second-derivative cross terms such as `Vρσ` (``\frac{∂^2e}{∂ρ∂σ}``).
"""
function kernel_terms end

#
# LDA
#
function potential_terms(func::Functional{:lda}, ρ::AbstractMatrix{T}) where {T}
    @assert has_energy(func)  # Otherwise custom implementation of this function needed
    s_ρ, n_p = size(ρ)
    TT = promote_type(T, parameter_type(func))

    e = mapreduce(hcat, ρ[:, i] for i = 1:n_p ) do ρi
        energy(func,ρi)
    end
    Vρ = mapreduce(hcat, ρ[:, i] for i = 1:n_p ) do ρi
        ForwardDiff.gradient(ρ -> energy(func,ρ), ρi)
    end
    (; e, Vρ)
end

function kernel_terms(func::Functional{:lda}, ρ::AbstractMatrix{T}) where {T}
    @assert has_energy(func)
    s_ρ, n_p = size(ρ)

    e = mapreduce(hcat, ρ[:, i] for i = 1:n_p ) do ρi
        energy(func,ρi)
    end
    Vρ = mapreduce(hcat, ρ[:, i] for i = 1:n_p ) do ρi
        ForwardDiff.gradient(ρ -> energy(func,ρ), ρi)
    end
    Vρρ = mapreduce(hcat, ρ[:, i] for i = 1:n_p ) do ρi
        ForwardDiff.hessian(ρ -> energy(func,ρ), ρi)
    end
    (; e, Vρ, Vρρ)
end

function energy(func::Functional{:lda}, ρ::AbstractVector{T}) where {T}
    length(ρ) == 1 || error("Multiple spins not yet implemented for fallback functionals")
    ρtotal = ρ[1]
    if ρtotal ≤ threshold_ρ(func, T)
        zero(T)
    else
        energy(func, ρtotal)
    end
end

#
# GGA
#
function potential_terms(func::Functional{:gga}, ρ::AbstractMatrix{T},
                         σ::AbstractMatrix{U}) where {T,U}
    @assert has_energy(func)  # Otherwise custom implementation of this function needed
    s_ρ, n_p = size(ρ)
    s_σ = size(σ, 1)

    e = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        energy(func, ρi, σi)
    end
    Vρ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        ForwardDiff.gradient(ρ -> energy(func, ρ, σi), ρi)
    end
    Vσ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        ForwardDiff.gradient(σ -> energy(func, ρi, σ), σi)
    end
    (; e, Vρ, Vσ)
end

function kernel_terms(func::Functional{:gga}, ρ::AbstractMatrix{T},
                      σ::AbstractMatrix{U}) where {T,U}
    @assert has_energy(func)  # Otherwise custom implementation of this function needed
    s_ρ, n_p = size(ρ)
    s_σ = size(σ, 1)
    TT = promote_type(T, U, parameter_type(func))

    e   = similar(ρ, TT, n_p)
    Vρ  = similar(ρ, TT, s_ρ, n_p)
    Vσ  = similar(ρ, TT, s_σ, n_p)
    Vρρ = similar(ρ, TT, s_ρ, s_ρ, n_p)
    Vρσ = similar(ρ, TT, s_ρ, s_σ, n_p)
    Vσσ = similar(ρ, TT, s_σ, s_σ, n_p)

    # TODO Needed to make forward-diff work with !isbits floating-point types (e.g. BigFloat)
    Vρ  .= zero(TT)
    Vσ  .= zero(TT)
    Vρρ .= zero(TT)
    Vρσ .= zero(TT)
    Vσσ .= zero(TT)

    e = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        energy(func, ρi, σi)
    end
    Vρ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        ForwardDiff.gradient(ρ -> energy(func, ρ, σi), ρi)
    end
    Vσ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        ForwardDiff.gradient(σ -> energy(func, ρi, σ), σi)
    end
    Vρρ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        ForwardDiff.hessian(ρ -> energy(func, ρ, σi), ρi)
    end
    Vρσ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        # mixed Hessian
        ForwardDiff.jacobian(
            σ -> ForwardDiff.gradient(ρ -> energy(func, ρ, σ), ρi), σi)
    end
    Vσσ = mapreduce(hcat, (ρ[:, i], σ[:, i]) for i = 1:n_p ) do (ρi, σi)
        ForwardDiff.hessian(σ -> energy(func, ρi, σ), σi)
    end

    (; e, Vρ, Vσ, Vρρ, Vρσ, Vσσ)
end

function energy(func::Functional{:gga}, ρ::AbstractVector{T},
                σ::AbstractVector{U}) where {T,U}
    length(ρ) == 1 || error("Multiple spins not yet implemented for fallback functionals")
    @assert length(ρ) == 1

    ρtotal = ρ[1]
    σtotal = σ[1]
    if ρtotal ≤ threshold_ρ(func, T)
        zero(promote_type(T, U, parameter_type(func)))
    else
        σstable = max(σtotal, threshold_σ(func, U))
        energy(func, ρtotal, σstable)
    end
end
