include("Misc.jl")
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

function BathGate(s1::Index, s2::Index, τ::Float64, ϵₛ::Float64, ϵ::Float64, γ::Float64, Γ::Float64, V::Float64, T::Float64, L::Int64)
    ρ = FermiDirac(ϵ, V, T)
    κ = √(Γ*γ/(2π))
    
    h = -im*γ*ρ*op("I",s1)*op("I",s2)
      + (ϵ-im*γ*(0.5-ρ))*op("nᴾ", s1)*op("I", s2)
      - (ϵ+im*γ*(0.5-ρ))*op("nᴬ", s1)*op("I", s2)
      + im*γ*ρ*op("b†ᴾ * b†ᴬ", s1)*op("I", s2)
      + im*γ*(1-ρ)*op("bᴾ * bᴬ", s1)*op("I", s2)
      + (ϵₛ/(2L))*op("I", s1)*(op("nᴾ", s2)-op("nᴬ", s2))
      + κ*op("b†ᴾ * JWᴬ", s1)*op("bᴾ", s2)
      + κ*op("bᴾ * JWᴬ", s1)*op("b†ᴾ", s2)
      - κ*op("bᴬ", s1)*op("JWᴾ * b†ᴬ", s2)
      - κ*op("b†ᴬ", s1)*op("JWᴾ * bᴬ", s2)
    return exp(-im*τ*h)
end

function SystemGate(s1::Index, s2::Index, τ::Float64, ϵₛ::Float64, t::Float64, U::Float64)
    h = 0.5*ϵₛ*(op("nᴾ", s1)-op("nᴬ", s1))*op("I", s2)
      + 0.5*ϵₛ*op("I", s1)*(op("nᴾ", s2)-op("nᴬ", s2))
      - t*op("b†ᴾ * JWᴬ", s1)*op("bᴾ", s2)
      - t*op("bᴾ * JWᴬ", s1)*op("b†ᴾ", s2)
      + t*op("b†ᴬ", s1)*op("JWᴾ * bᴬ", s2)
      + t*op("bᴬ", s1)*op("JWᴾ * b†ᴬ", s2)
      + U*op("nᴾ", s1)*op("nᴾ", s2)
      - U*op("nᴬ", s1)*op("nᴬ", s2)
    return exp(-im*τ*h)
end

function TimeEvolutionOperator(sites::Vector{<:Index}, δτ::Float64, ϵ₀::Float64, t::Float64, U::Float64, ϵ::Array{Float64}, γ::Array{Float64}, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, L::Int64, D::Int64)

    LeftBathGates = ITensor[]
    for n in 1:L
        s1 = sites[L-n+1]
        s2 = sites[L-n+2]
        push!(LeftBathGates, BathGate(s1, s2, δτ/2, ϵ₀, ϵ[n], γ[n], Γ, -ΔV/2, Tₗ, L))
        push!(LeftBathGates, op("SWAP", s1, s2))
    end
    append!(LeftBathGates, reverse(LeftBathGates))

    RightBathGates = ITensor[]
    for n in 1:L
        s1 = sites[L+D+n-1]
        s2 = sites[L+D+n]
        push!(RightBathGates, BathGate(s2, s1, δτ/2, ϵ₀, ϵ[n], γ[n], Γ, ΔV/2, Tᵣ, L))
        push!(RightBathGates, op("SWAP", s1, s2))
    end
    append!(RightBathGates, reverse(RightBathGates))

    SystemGates = ITensor[]
    for n in 1:D-1
        s1 = sites[L+n]
        s2 = sites[L+n+1]
        push!(SystemGates, SystemGate(s1, s2, δτ/2, ϵ₀, t, U))
    end

    TimeEvolution = ITensor[]
    append!(TimeEvolution, LeftBathGates)
    append!(TimeEvolution, SystemGates)
    append!(TimeEvolution, RightBathGates)
    append!(TimeEvolution, reverse(SystemGates))

    return TimeEvolution
end

#Machine parameters
const ϵ₀ = 1/8 #1.
const t = 0. #1.
const U = 0. #1.2
const Γ = 1/2 #6.
const ΔV = 1/8 #1.
const Tₗ = 1/8 #10.
const Tᵣ = 1/8 #1.
const W = 1. #8.
const W° = 1/2 #4.
const L₁ = 4 
const L₂ = 2 
const L = L₁+L₂ 
const D = 1

const δτ = 0.01
const N = 1
const n_print = 1
const χ = 40

let
    println("Expected Values: Jₚ=0.0062268975080337265, Jₕ=0.0006487902180691001\n")
    ITensors.enable_debug_checks()

    sites = siteinds("SuperFermion", 2*L+D)
    I_vacc = MPS(sites, ["Vacuum" for n in 1:length(sites)])
    ρ̂ = MPS(sites, ["NormalizedVacuum" for n in 1:length(sites)])

    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    Ĵₚ = ParticleCurrentOperator(sites, ϵ, γ, -ΔV/2, Tₗ)
    Ĵₕ = EnergyCurrentOperator(sites, ϵ, γ, Γ, -ΔV/2, Tₗ)
    println("t,Jₚ,Jₕ")
    println(0., ",", real(inner(I_vacc', Ĵₚ, ρ̂)), ",", real(inner(I_vacc', Ĵₕ, ρ̂)),",1")
    
    Û = TimeEvolutionOperator(sites, δτ, ϵ₀, t, U, ϵ, γ, Γ, ΔV, Tₗ, Tᵣ, L, D)
    @show Û
    #=
    for n in 1:N
        ρ̂ = apply(Û, ρ̂; maxdim=χ)
        trace = real(dot(I_vacc, ρ̂))
        if n%n_print == 0
            println(round(n*δτ, digits=2),",", real(inner(I_vacc', Ĵₚ, ρ̂))/trace,",", real(inner(I_vacc', Ĵₕ, ρ̂))/trace,",",trace)
        else
            println(round(n*δτ, digits=2),",,,")
        end
    end
    =#

    return
end
