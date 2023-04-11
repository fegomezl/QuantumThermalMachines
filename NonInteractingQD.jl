using LinearAlgebra
using Roots

using DataFrames
using CSV

using LazyArrays
using FillArrays
using StaticArrays

function FermiDirac(ϵ::Float64, μ::Float64, T::Float64)
    return 1/(1+exp((ϵ-μ)/T))
end

function ExponentialRate(L::Int64, C::Float64)
    f(x) = x*(1.0-x^L)/(1.0-x+1e-12)-C
    return find_zero(f, (1.0, C), A42())
end

function BathSpectra(W::Float64, W°::Float64, L₁::Int64, L₂::Int64)
    Δϵ = 2.0*W°/(L₁-1)
    Φ = ExponentialRate(L₂÷2, (W-W°)/Δϵ)
    
    Range1 = 1:L₁÷2
    Range2 = 1:L₂÷2
    exponent = @. Φ^Range2
    aux = (W-W°)/(1-Φ^(L₂÷2))
    
    ϵ₁ = @. Δϵ*(Range1 - 0.5)
    ϵ₂ = @. W° + aux *(1.0 - exponent)
    
    ϵ = ApplyArray(vcat, reverse(-ϵ₂), reverse(-ϵ₁), ϵ₁, ϵ₂)

    γ₁ = Fill(Δϵ, L₁÷2)
    γ₂ = @. Δϵ*Φ^ϵ₂

    γ = ApplyArray(vcat, reverse(γ₂), reverse(γ₁), γ₁, γ₂)
    
    return ϵ, γ
end

function SystemBathHamiltonian(ϵ₀::Float64, Γ::Float64, ϵ, γ)
    Hₛ = [ϵ₀]
    Hₗ = Diagonal(ApplyArray(vcat, ϵ, ϵ))
    
    aux = sqrt(Γ/(2π))
    κ =  @. aux*sqrt(γ)
    Hₛₗ = ApplyArray(vcat, κ, κ)

    H₁ = ApplyArray(hcat, Hₗ, Hₛₗ)
    H₂ = ApplyArray(hcat, Hₛₗ', Hₛ)

    return ApplyArray(vcat, H₁, H₂)
end

function TunnelingRates(ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ, γ)
    ρₗ = @. γ*FermiDirac(ϵ, ΔV/2, Tₗ)
    ρᵣ = @. γ*FermiDirac.(ϵ, -ΔV/2, Tᵣ)
    
    Γ₊ = ApplyArray(vcat, ρₗ, ρᵣ, zeros(SVector{1}))
    Γ₋ = ApplyArray(vcat, γ.-ρₗ, γ.-ρᵣ, zeros(SVector{1}))
    
    return Diagonal(Γ₊), Diagonal(Γ₋)
end

function Liouvillian(ϵ₀::Float64, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ, γ)
    H = SystemBathHamiltonian(ϵ₀, Γ, ϵ, γ)
    Γ₊, Γ₋ = TunnelingRates(ΔV, Tₗ, Tᵣ, ϵ, γ)
    Ω = (Γ₋ - Γ₊)/2
    
    A1 = @. H-im*Ω
    B1 = @. im*Γ₊
    L₁ = ApplyArray(hcat, A1, B1)
    
    A2 = @. im*Γ₋
    B2 = @. H+im*Ω
    L₂ = ApplyArray(hcat, A2,  B2)
    
    return ApplyArray(vcat, L₁, L₂)
end

function RunMachine(ϵ₀::Float64, W::Float64, W°::Float64, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, L₁::Int64, L₂::Int64)
    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    L = Liouvillian(ϵ₀, Γ, ΔV, Tₗ, Tᵣ, ϵ, γ) 
    
    λ, V = eigen(materialize(L),permute=false,scale=false)
    V⁻¹ = inv(materialize(V))
    D = @. 0.5*(sign(imag(λ))+1)
    Corr = V*Diagonal(D)*V⁻¹

    aₖaₖ= diag(Corr)[1:L₁+L₂]
    aₖc = Corr[1:L₁+L₂, 1+2(L₁+L₂)]
    caₖ = Corr[1+2(L₁+L₂), 1:L₁+L₂]

    A = @. FermiDirac(ϵ, ΔV/2, Tₗ) - real(aₖaₖ)
    B = @. real(aₖc) + real(caₖ)

    Jₚ = @. sum(γ*A)
    Jₕ = @. sum(γ*ϵ*A - √(Γ/(8π))*(γ^1.5)*B)

    return Jₚ, Jₕ
end

#Machine parameters
const ϵ₀ = 1/8
const ΔV = 1/8
const Tₗ = 1/8
const Tᵣ = 1/8
const W = 1.0
const W° = 1/2
const L₁ = 160
const L₂ = 40

#Sweep in tunneling strength
Γ = LinRange(0.0, 0.5, 20)
J = RunMachine.(ϵ₀, W, W°, Γ, ΔV, Tₗ, Tᵣ, L₁, L₂)
J = reinterpret(reshape, Float64, J)
data = DataFrame(Γ=Γ, Jₚ=J[1,:], Jₕ=J[2,:])
CSV.write("data.csv", data)