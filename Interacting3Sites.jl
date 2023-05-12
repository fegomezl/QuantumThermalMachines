using LinearAlgebra
using Roots
using ITensors

# Fermi-Dirac distribution
function FermiDirac(ϵ::Float64, μ::Float64, T::Float64)
    return 1/(1+exp((ϵ-μ)/T))
end

# Find the spacing in the logaritmic region:
# Solve C from ϵ(n)=A+B*Cⁿ
# With ϵ(0) = W°-Δϵ
#      ϵ(1) = W°
#      ϵ(L+1) = W
# RHS = (W-W°)/Δϵ
function ExponentialRate(L::Int64, RHS::Float64)
    f(x) = x*(1-x^L)/(1-x+1e-9)-RHS
    return find_zero(f, (1, RHS), Bisection()) 
    #RHS is a high bound for the solution and allows to use bisection
end

# Get energies and spacings for a uniform bath with 
# energy window [-W,W] and enhanced resolution in
# [-W°,W°] with L₁ points inside and L₂ outside
function BathSpectra(W::Float64, W°::Float64, L₁::Int64, L₂::Int64)

    # Uniform energies
    Δϵ = 2*W°/(L₁-1) 
    ϵ₁ = [1:1:L₁÷2;]
    ϵ₁ = Δϵ.*(ϵ₁.-0.5)

    # Logarithmic energies
    Φ = ExponentialRate(L₂÷2, (W-W°)/Δϵ)
    ϵ₂ = [1:1:L₂÷2;]
    ϵ₂ = W° .+ (W-W°).*(1.0.-Φ.^ϵ₂)./(1-Φ^(L₂÷2))

    ϵ = vcat(reverse(-ϵ₂), reverse(-ϵ₁), ϵ₁, ϵ₂)

    # Uniform spacing
    γ₁ = Δϵ.*ones(L₁÷2)

    # Logarithmic spacing
    γ₂ = [1:1:L₂÷2;]
    γ₂ = Δϵ.*Φ.^ϵ₂

    γ = vcat(reverse(γ₂), reverse(γ₁), γ₁, γ₂)

    return ϵ, γ
end

include("SuperFermionSite.jl")

function ParticleCurrentOperator(sites::Vector{<:Index}, ϵ::Array{Float64}, γ::Array{Float64}, V::Float64, T::Float64)
    Ĵₚ = OpSum()
    for n in eachindex(ϵ)
        Ĵₚ += γ[n]*FermiDirac(ϵ[n],V,T),"I",n
        Ĵₚ += -γ[n],"nᴾ",n
    end
    return MPO(Ĵₚ, sites)
end

function EnergyCurrentOperator(sites::Vector{<:Index}, ϵ::Array{Float64}, γ::Array{Float64}, Γ::Float64, V::Float64, T::Float64)
    Ĵₕ = OpSum()
    SystemIndex = length(ϵ)+1
    for n in eachindex(ϵ)
        Ĵₕ += γ[n]*ϵ[n]*FermiDirac(ϵ[n],V,T),"I",n
        Ĵₕ += -γ[n]*ϵ[n],"nᴾ",n
        JW = ()
        for m in n+1:SystemIndex-1
            JW = (JW..., "JWᴾ * JWᴬ",m)
        end
        tₙ = √(Γ/(8π))*γ[n]^1.5
        Ĵₕ += -tₙ,"b†ᴾ * JWᴬ",n,JW...,"bᴾ",SystemIndex
        Ĵₕ += -tₙ,"b†ᴾ",SystemIndex,JW...,"JWᴬ * bᴾ",n
    end
    return MPO(Ĵₕ, sites)
end

#Machine parameters
const ϵ₀ = 1.
const t = 1.
const U = 0. #1.2
const Γ = 6.
const ΔV = 1.
const Tₗ = 10.
const Tᵣ = 1.
const W = 8.
const W° = 4.
const L₁ = 40
const L₂ = 10

N = L₁+L₂+1
sites = siteinds("SuperFermion", N)

I_vacc = MPS(sites, ["Vacuum" for n in 1:N])
ρ₀ = MPS(sites, ["NormalizedVacuum" for n in 1:N])

ϵ, γ = BathSpectra(W, W°, L₁, L₂)
Ĵₚ = ParticleCurrentOperator(sites, ϵ, γ, -ΔV/2, Tₗ)
Ĵₕ = EnergyCurrentOperator(sites, ϵ, γ, Γ, -ΔV/2, Tₗ)

println(inner(I_vacc', Ĵₚ, ρ₀))
println(sum(γ.*(FermiDirac.(ϵ,-ΔV/2,Tₗ).-0.5)))

println(inner(I_vacc', Ĵₕ, ρ₀))
println(sum(γ.*ϵ.*(FermiDirac.(ϵ,-ΔV/2,Tₗ).-0.5)))
