include("Misc.jl")
include("SuperFermionSite.jl")

function BathParticleCurrentOperator(sites::Vector{<:Index}, ϵ::Array{Float64}, γ::Array{Float64}, n0::Int64, V::Float64, T::Float64)
    Ĵₚ = OpSum()
    for n in eachindex(ϵ)
        Ĵₚ += γ[n]*FermiDirac(ϵ[n],V,T),"I",n0+n
        Ĵₚ += -γ[n],"nᴾ",n0+n
    end
    return MPO(Ĵₚ, sites)
end

function SystemParticleCurrentOperator(sites::Vector{<:Index}, index1::Int64, index2::Int64, t::Float64)
    Ĵₚ = OpSum()
    Ĵₚ += im*t,"b†ᴾ * JWᴬ",index1,"bᴾ",index2
    Ĵₚ += -im*t,"b†ᴾ",index2,"JWᴬ * bᴾ",index2
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
    
    h = -im*γ*ρ*op("I",s1)*op("I",s2) +
        (ϵ-im*γ*(0.5-ρ))*op("nᴾ", s1)*op("I", s2) -
        (ϵ+im*γ*(0.5-ρ))*op("nᴬ", s1)*op("I", s2) +
        im*γ*ρ*op("b†ᴾ * b†ᴬ", s1)*op("I", s2) +
        im*γ*(1-ρ)*op("bᴾ * bᴬ", s1)*op("I", s2) +
        (ϵₛ/L)*op("I", s1)*(op("nᴾ", s2)-op("nᴬ", s2)) +
        κ*op("b†ᴾ * JWᴬ", s1)*op("bᴾ", s2) +
        κ*op("bᴾ * JWᴬ", s1)*op("b†ᴾ", s2) -
        κ*op("bᴬ", s1)*op("JWᴾ * b†ᴬ", s2) -
        κ*op("b†ᴬ", s1)*op("JWᴾ * bᴬ", s2)
    return exp(-im*τ*h)
end

function SystemGate(s1::Index, s2::Index, τ::Float64, ϵₛ::Float64, t::Float64, U::Float64)
    h = 0.5*ϵₛ*(op("nᴾ", s1)-op("nᴬ", s1))*op("I", s2) + 
        0.5*ϵₛ*op("I", s1)*(op("nᴾ", s2)-op("nᴬ", s2)) - 
        t*op("b†ᴾ * JWᴬ", s1)*op("bᴾ", s2) - 
        t*op("b†ᴾ", s2)*op("JWᴬ * bᴾ", s1) + 
        t*op("b†ᴬ", s1)*op("JWᴾ * bᴬ", s2) + 
        t*op("b†ᴬ * JWᴾ", s2)*op("bᴬ", s1) + 
        U*op("nᴾ", s1)*op("nᴾ", s2) - 
        U*op("nᴬ", s1)*op("nᴬ", s2)
    return exp(-im*τ*h)
end

function TimeEvolutionOperator(sites::Vector{<:Index}, δτ::Float64, ϵ₀::Float64, t::Float64, U::Float64, ϵ::Array{Float64}, γ::Array{Float64}, Γ::Float64, Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, L::Int64, D::Int64)

    LeftBathGates = ITensor[]
    for n in 1:L
        s1 = sites[L-n+1]
        s2 = sites[L-n+2]
        push!(LeftBathGates, BathGate(s1, s2, δτ/2, ϵ₀, ϵ[n], γ[n], Γ, Vₗ, Tₗ, L))
        push!(LeftBathGates, op("SWAP", s1, s2))
    end
    for n in 1:L
        s1 = sites[n]
        s2 = sites[n+1]
        push!(LeftBathGates, BathGate(s2, s1, δτ/2, ϵ₀, ϵ[L-n+1], γ[L-n+1], Γ, Vₗ, Tₗ, L))
        push!(LeftBathGates, op("SWAP", s2, s1))
    end

    RightBathGates = ITensor[]
    for n in 1:L
        s1 = sites[L+D+n-1]
        s2 = sites[L+D+n]
        push!(RightBathGates, BathGate(s2, s1, δτ/2, ϵ₀, ϵ[n], γ[n], Γ, Vᵣ, Tᵣ, L))
        push!(RightBathGates, op("SWAP", s2, s1))
    end
    for n in 1:L
        s1 = sites[2L+D-n]
        s2 = sites[2L+D-n+1]
        push!(RightBathGates, BathGate(s1, s2, δτ/2, ϵ₀, ϵ[L-n+1], γ[L-n+1], Γ, Vᵣ, Tᵣ, L))
        push!(RightBathGates, op("SWAP", s1, s2))
    end

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
const ϵ₀ = 1.
const t = 1.
const U = 1.2 
const W = 8.
const W° = 4.
const Γ = 3.
const Vₗ = -0.5
const Vᵣ = 0.5
const Tₗ = 10.
const Tᵣ = 1.
const L₁ = 4 
const L₂ = 2 
const L = L₁+L₂ 
const D = 3

const δτ = 0.05
const N = 1000
const n_print = 50
const χ = 40

let
    sites = siteinds("SuperFermion", 2*L+D)
    I_vacc = MPS(sites, ["Vacuum" for n in 1:length(sites)])
    ρ̂ = MPS(sites, ["NormalizedVacuum" for n in 1:length(sites)])

    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    Ĵₚᴸ = BathParticleCurrentOperator(sites, ϵ, γ, 0, Vₗ, Tₗ)
    Ĵₚᴿ = BathParticleCurrentOperator(sites, ϵ, γ, L+D, Vᵣ, Tᵣ)
    Ĵₚ¹² = SystemParticleCurrentOperator(sites, L+1, L+2, t)
    Ĵₚ²³ = SystemParticleCurrentOperator(sites, L+2, L+3, t)
    Û = TimeEvolutionOperator(sites, δτ, ϵ₀, t, U, ϵ, γ, Γ, Vₗ, Vᵣ, Tₗ, Tᵣ, L, D)

    println("t,Jₚᴸ,Jₚᴿ,Jₚ¹²,Jₚ²³, trace")
    println(0., ",", 
            real(inner(I_vacc', Ĵₚᴸ, ρ̂)), ",", 
            real(inner(I_vacc', Ĵₚᴿ, ρ̂)), ",", 
            real(inner(I_vacc', Ĵₚ¹², ρ̂)), ",", 
            real(inner(I_vacc', Ĵₚ²³, ρ̂)), ",", 
            "1") 
    for n in 1:N
        ρ̂ = apply(Û, ρ̂; maxdim=χ)
        if n%n_print == 0
            trace = real(dot(I_vacc, ρ̂))
            println(round(n*δτ, digits=2), ",", 
                    real(inner(I_vacc', Ĵₚᴸ, ρ̂))/trace, ",", 
                    real(inner(I_vacc', Ĵₚᴿ, ρ̂))/trace, ",", 
                    real(inner(I_vacc', Ĵₚ¹², ρ̂))/trace, ",", 
                    real(inner(I_vacc', Ĵₚ²³, ρ̂))/trace, ",", 
                    trace) 
        else
            println(round(n*δτ, digits=2),",,,")
        end
    end
    return
end
