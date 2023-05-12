using LinearAlgebra
using Roots

using DataFrames
using CSV

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

# Calculate the system-lead hamiltonian:
# H = [ Hₗ   Hₛₗ ]
#     [ Hₛₗ' Hₛ  ]
# Hₗ = diag(ϵ₁,ϵ₂,ϵ₃,ϵ₄...) ϵᵢ: Energies from the leads
# Hₛ = [ϵ₀] ~ Single site system
# Hₛₗ = [κ₁,κ₂,κ₃,κ₄...] κᵢ = √(Γγᵢ/2π) γᵢ: Spacings from the leads
function SystemLeadHamiltonian(ϵ₀::Float64, Γ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    Hₛ = [ϵ₀]
    Hₗ = Diagonal(vcat(ϵ, ϵ))

    κ = sqrt.(Γ.*γ./(2π))
    Hₛₗ = vcat(κ, κ)

    H₁ = hcat(Hₗ, Hₛₗ)
    H₂ = hcat(Hₛₗ', Hₛ)
    return vcat(H₁, H₂)
end

# Calculate the tunneling rates 
# Γ₊= diag(γ₁ρ₁,γ₂ρ₂,γ₃ρ₃...0) ρᵢ: Fermi-Dirac distribution on given lead
# Γ₋= diag(γ₁(1-ρ₁),γ₂(1-ρ₂),γ₃(1-ρ₃)...0) ρᵢ: Fermi-Dirac distribution on given lead
# The zeroes correspond to system sites
function TunnelingRates(ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    ρₗ = FermiDirac.(ϵ, ΔV/2, Tₗ)
    ρᵣ = FermiDirac.(ϵ, -ΔV/2, Tᵣ)
    Γ₊ = vcat(γ.*ρₗ, γ.*ρᵣ, [0.0])
    Γ₋ = vcat(γ.*(1.0.-ρₗ), γ.*(1.0.-ρᵣ), [0.0])
    return Diagonal(Γ₊), Diagonal(Γ₋)
end

# Create the liouvillian for the whole system:
# L = [H-iΩ  iΓ₊ ]
#     [iΓ₋   H+iΩ]
function Liouvillian(ϵ₀::Float64, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    H = SystemLeadHamiltonian(ϵ₀, Γ, ϵ, γ)
    Γ₊, Γ₋ = TunnelingRates(ΔV, Tₗ, Tᵣ, ϵ, γ)
    Ω = (Γ₋ - Γ₊)/2

    L₁ = hcat(H-im*Ω, im*Γ₊ )
    L₂ = hcat(im*Γ₋,  H+im*Ω)
    return vcat(L₁, L₂)
end

# Perform the simulation for a given set of parameters and return the corresponding
# particle and energy currents Jₚ, Jₕ.
# Parameters:
#   -ϵ₀: Single site energy of the system
#   -W: Energy window for the bath
#   -W°: Enhanced resolution energy window for the bath
#   -Γ: Tunneling rate for the bath
#   -ΔV: Potential difference
#   -Tₗ: Left temperature
#   -Tᵣ: Right temperature
#   -L₁: Points inside the enhanced resolution energy window.
#   -L₂: Points outside the enhanced resolution energy window.
function RunMachine(ϵ₀::Float64, W::Float64, W°::Float64, Γ::Float64, ΔV::Float64, Tₗ::Float64, Tᵣ::Float64, L₁::Int64, L₂::Int64)
    print("Running for Γ="*string(Γ)*" ... ")

    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    L = Liouvillian(ϵ₀, Γ, ΔV, Tₗ, Tᵣ, ϵ, γ) 

    λ, V = eigen(L)
    V⁻¹ = inv(V)

    #Dᵢ = { 1 ; imag(λᵢ)>0
    #     { 0 ; imag(λᵢ)<0
    D = (sign.(imag.(λ)).+1)./2
    Corr = V*Diagonal(D)*V⁻¹

    aₖaₖ= diag(Corr)[1:L₁+L₂]
    aₖc = Corr[1:L₁+L₂, 1+2(L₁+L₂)]
    caₖ = Corr[1+2(L₁+L₂), 1:L₁+L₂]

    A = FermiDirac.(ϵ, ΔV/2, Tₗ).-real.(aₖaₖ)
    B = real.(aₖc)+real.(caₖ)

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
