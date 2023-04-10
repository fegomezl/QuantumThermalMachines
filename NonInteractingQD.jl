using LinearAlgebra
using Roots

using DataFrames
using CSV

function FermiDirac(ϵ::Float64, μ::Float64, T::Float64)
    return 1/(1+exp((ϵ-μ)/T))
end

function ExponentialRate(L::Int64, C::Float64)
    f(x) = x*(1-x^L)/(1-x+1e-9)-C
    return find_zero(f, (1, C), Bisection())
end

function BathSpectra(W::Float64, W°::Float64, L₁::Int64, L₂::Int64)
    Δϵ = 2*W°/(L₁-1)
    Φ = ExponentialRate(L₂÷2, (W-W°)/Δϵ)

    ϵ₁ = [1:1:L₁÷2;]
    ϵ₁ = Δϵ.*(ϵ₁.-0.5)

    ϵ₂ = [1:1:L₂÷2;]
    ϵ₂ = W° .+ (W-W°).*(1.0.-Φ.^ϵ₂)./(1-Φ^(L₂÷2))

    ϵ = vcat(reverse(-ϵ₂), reverse(-ϵ₁), ϵ₁, ϵ₂)

    γ₁ = Δϵ.*ones(L₁÷2)

    γ₂ = [1:1:L₂÷2;]
    γ₂ = Δϵ.*Φ.^ϵ₂

    γ = vcat(reverse(γ₂), reverse(γ₁), γ₁, γ₂)

    return ϵ, γ
end

function SystemBathHamiltonian(ϵ₀::Float64, Γ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    Hₛ = [ϵ₀]
    Hₗ = Diagonal(vcat(ϵ, ϵ))

    κ = sqrt.(Γ.*γ./(2π))
    Hₛₗ = vcat(κ, κ)

    H₁ = hcat(Hₗ, Hₛₗ)
    H₂ = hcat(Hₛₗ', Hₛ)
    return vcat(H₁, H₂)
end


function TunnelingRates(ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    ρₗ = FermiDirac.(ϵ, ΔV, Tₗ)
    ρᵣ = FermiDirac.(ϵ, 0.0, Tᵣ)
    Γ₊ = vcat(γ.*ρₗ, γ.*ρᵣ, [0.0])
    Γ₋ = vcat(γ.*(1.0.-ρₗ), γ.*(1.0.-ρᵣ), [0.0])
    return Diagonal(Γ₊), Diagonal(Γ₋)
end

function Liouvillian(ϵ₀::Float64, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    H = SystemBathHamiltonian(ϵ₀, Γ, ϵ, γ)
    Γ₊, Γ₋ = TunnelingRates(ΔV, Tₗ, Tᵣ, ϵ, γ)
    Ω = (Γ₋ - Γ₊)/2

    L₁ = hcat(H-im*Ω, im*Γ₊ )
    L₂ = hcat(im*Γ₋,  H+im*Ω)
    return vcat(L₁, L₂)
end

function RunMachine(ϵ₀::Float64, W::Float64, W°::Float64, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, L₁::Int64, L₂::Int64)
    print("Running for Γ="*string(Γ)*" ... ")

    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    L = Liouvillian(ϵ₀, Γ, ΔV, Tₗ, Tᵣ, ϵ, γ) 

    λ, V = eigen(L)
    V⁻¹ = inv(V)
    D = (sign.(imag.(λ)).+1)./2
    Corr = V*Diagonal(D)*V⁻¹

    aₖaₖ= diag(Corr)[1:L₁+L₂]
    aₖc = Corr[1:L₁+L₂, 1+2(L₁+L₂)]
    caₖ = Corr[1+2(L₁+L₂), 1:L₁+L₂]

    A = FermiDirac.(ϵ, ΔV, Tₗ).-real.(aₖaₖ)
    B = (real.(aₖc)+real.(caₖ))./2

    Jₚ = sum(γ.*A)
    Jₕ = sum(γ.*ϵ.*A .- √(Γ/(8π)).*(γ.^1.5).*B)

    println("Done!")
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
