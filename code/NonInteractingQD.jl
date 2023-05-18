#!/usr/bin/julia
include("Misc.jl")

# Calculate the system-lead hamiltonian:
# H = [ Hₗ   Hₛₗ ]
#     [ Hₛₗ' Hₛ  ]
# Hₗ = diag(ϵ₁,ϵ₂,ϵ₃,ϵ₄...) ϵᵢ: Energies from the leads
# Hₛ = [ϵ₀] ~ Single site system
# Hₛₗ = [κ₁,κ₂,κ₃,κ₄...] κᵢ = √(Γγᵢ/2π) γᵢ: Spacings from the leads
function SystemLeadHamiltonian(ϵ₀::Float64, Γ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    Hₛ = [ϵ₀]
    Hₗ = Diagonal(vcat(ϵ, ϵ))

    κ = @. √(Γ/(2π))*√γ
    Hₛₗ = vcat(κ, κ)

    H₁ = hcat(Hₗ, Hₛₗ)
    H₂ = hcat(Hₛₗ', Hₛ)
    return vcat(H₁, H₂)
end

# Calculate the tunneling rates 
# Γ₊= diag(γ₁ρ₁,γ₂ρ₂,γ₃ρ₃...0) ρᵢ: Fermi-Dirac distribution on given lead
# Γ₋= diag(γ₁(1-ρ₁),γ₂(1-ρ₂),γ₃(1-ρ₃)...0) ρᵢ: Fermi-Dirac distribution on given lead
# The zeroes correspond to system sites
function TunnelingRates(Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    fₗ = FermiDirac.(ϵ, Vₗ, Tₗ)
    fᵣ = FermiDirac.(ϵ, Vᵣ, Tᵣ)
    Γ₊ = vcat(γ.*fₗ, γ.*fᵣ, [0.0])
    Γ₋ = vcat(γ.*(1.0.-fₗ), γ.*(1.0.-fᵣ), [0.0])
    return Diagonal(Γ₊), Diagonal(Γ₋)
end

# Create the liouvillian for the whole system:
# L = [H-iΩ  iΓ₊ ]
#     [iΓ₋   H+iΩ]
function Liouvillian(ϵ₀::Float64, Γ::Float64, Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    H = SystemLeadHamiltonian(ϵ₀, Γ, ϵ, γ)
    Γ₊, Γ₋ = TunnelingRates(Vₗ, Vᵣ, Tₗ, Tᵣ, ϵ, γ)
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
function RunMachine(ϵ₀::Float64, W::Float64, W°::Float64, Γ::Float64, Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, L₁::Int64, L₂::Int64)

    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    L = Liouvillian(ϵ₀, Γ, Vₗ, Vᵣ, Tₗ, Tᵣ, ϵ, γ) 

    # Diagonalize L = V⁻¹λV
    λ, V = eigen(L,permute=false,scale=false)
    V⁻¹ = inv(V)

    #Dᵢ = { 1 ; imag(λᵢ)>0
    #     { 0 ; imag(λᵢ)<0
    D = @. (sign(imag(λ))+1)/2
    Corr = V*Diagonal(D)*V⁻¹

    #Calculate the expected values for the currents
    aₖaₖ= diag(Corr)[1:L₁+L₂]
    aₖc = Corr[1:L₁+L₂, 1+2(L₁+L₂)]
    caₖ = Corr[1+2(L₁+L₂), 1:L₁+L₂]

    A = @. FermiDirac(ϵ, ΔV/2, Tₗ)-real(aₖaₖ)
    B = @. real(aₖc)+real(caₖ)

    Jₚ = sum(@. γ*A)
    Jₕ = sum(@. γ*ϵ*A - √(Γ/(8π))*(γ^1.5)*B)

    return Jₚ, Jₕ
end

#Machine parameters
const ϵ₀ = 1/8
const Vₗ = 1/16
const Vᵣ = -1/16
const ΔV = 1/8
const Tₗ = 1/8
const Tᵣ = 1/8
const W = 1.0
const W° = 1/2
const Γ₀ = 1/8

#=Sweep in tunneling strength
println("Sweep in tunneling strength")
for L in [50, 100, 200]
    println("L=$L")
    L₁ = Int(0.8*L)
    L₂ = Int(0.2*L)
    Γ = LinRange(0.0, 0.5, 20)
    J = RunMachine.(ϵ₀, W, W°, Γ, Vₗ, Vᵣ, Tₗ, Tᵣ, L₁, L₂)
    J = reinterpret(reshape, Float64, J)    
    data = DataFrame(Γ=Γ, Jₚ=J[1,:], Jₕ=J[2,:])
    CSV.write("results/gamma_"*string(L)*".csv", data)
end
=#

#=Currents vs Temperature
println("Currents vs Temperature")
for L in [50, 100, 200]
    println("L=$L")
    L₁ = Int(0.8*L)
    L₂ = Int(0.2*L)
    T = exp10.(range(-4.0,1.0,100))
    J = RunMachine.(ϵ₀, W, W°, Γ₀, Vₗ, Vᵣ, T, T, L₁, L₂)
    J = reinterpret(reshape, Float64, J)
    data = DataFrame(T=T, Jₚ=J[1,:], Jₕ=J[2,:])
    CSV.write("results/temp_"*string(L)*".csv", data)
end
=#

#=Currents vs single level energy
println("Currents vs single level energy")
for L in [20, 50, 100]
    println("L=$L")
    L₁ = Int(0.8*L)
    L₂ = Int(0.2*L)
    ϵ = LinRange(-1.0,1.0,100)
    J = RunMachine.(ϵ, W, W°, Γ₀, Vₗ, Vᵣ, Tₗ, Tᵣ, L₁, L₂)
    J = reinterpret(reshape, Float64, J)
    data = DataFrame(ϵ=ϵ, Jₚ=J[1,:], Jₕ=J[2,:])
    CSV.write("results/energy_"*string(L)*".csv", data)
end
=#
