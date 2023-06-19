#!/usr/bin/julia
include("Misc.jl")

"""
 Calculate the system-lead hamiltonian:
  H = [ Hₗ   Hₛₗ ]
      [ Hₛₗ' Hₛ  ]
 Hₗ = diag(ϵ₁,ϵ₂,ϵ₃,ϵ₄...) ϵᵢ: Energies from the leads
 Hₛ = [ϵ₀] ~ Single site system
 Hₛₗ = [κ₁,κ₂,κ₃,κ₄...] κᵢ = √(Γγᵢ/2π) γᵢ: Spacings from the leads
"""
function SystemLeadHamiltonian(ϵ₀::Float64, Γ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    Hₛ = [ϵ₀]
    Hₗ = Diagonal(vcat(ϵ, ϵ))

    κ = @. √(Γ/(2π))*√γ
    Hₛₗ = vcat(κ, κ)

    H₁ = hcat(Hₗ, Hₛₗ)
    H₂ = hcat(Hₛₗ', Hₛ)
    return vcat(H₁, H₂)
end

"""
 Calculate the tunneling rates 
 Γ₊= diag(γ₁ρ₁,γ₂ρ₂,γ₃ρ₃...0) ρᵢ: Fermi-Dirac distribution on given lead
 Γ₋= diag(γ₁(1-ρ₁),γ₂(1-ρ₂),γ₃(1-ρ₃)...0) ρᵢ: Fermi-Dirac distribution on given lead
 The zeroes correspond to system sites
"""
function TunnelingRates(Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    fₗ = FermiDirac.(ϵ, Vₗ, Tₗ)
    fᵣ = FermiDirac.(ϵ, Vᵣ, Tᵣ)
    Γ₊ = vcat(γ.*fₗ, γ.*fᵣ, [0.0])
    Γ₋ = vcat(γ.*(1.0.-fₗ), γ.*(1.0.-fᵣ), [0.0])
    return Diagonal(Γ₊), Diagonal(Γ₋)
end

"""
 Create the liouvillian for the whole system:
 L = [H-iΩ  iΓ₊ ]
     [iΓ₋   H+iΩ]
"""
function Liouvillian(ϵ₀::Float64, Γ::Float64, Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, ϵ::Array{Float64}, γ::Array{Float64})
    H = SystemLeadHamiltonian(ϵ₀, Γ, ϵ, γ)
    Γ₊, Γ₋ = TunnelingRates(Vₗ, Vᵣ, Tₗ, Tᵣ, ϵ, γ)
    Ω = (Γ₋ - Γ₊)/2

    L₁ = hcat(H-im*Ω, im*Γ₊ )
    L₂ = hcat(im*Γ₋,  H+im*Ω)
    return vcat(L₁, L₂)
end

"""
 Perform the simulation for a given set of parameters and return the corresponding
 particle and energy currents Jₚ, Jₕ.
 Parameters:
   -ϵ₀: Single site energy of the system
   -W: Energy window for the bath
   -W°: Enhanced resolution energy window for the bath
   -Γ: Tunneling rate for the bath
   -ΔV: Potential difference
   -Tₗ: Left temperature
   -Tᵣ: Right temperature
   -L₁: Points inside the enhanced resolution energy window.
   -L₂: Points outside the enhanced resolution energy window.
"""
function RunMachine(ϵ₀::Float64, W::Float64, W°::Float64, Γ::Float64, Vₗ::Float64, Vᵣ::Float64, Tₗ::Float64, Tᵣ::Float64, L₁::Int64, L₂::Int64)

    ϵ, γ = BathSpectra(W, W°, L₁, L₂)
    L = Liouvillian(ϵ₀, Γ, Vₗ, Vᵣ, Tₗ, Tᵣ, ϵ, γ) 

    # Diagonalize L = V⁻¹λV
    λ, V = eigen(L)
    V⁻¹ = inv(V)

    #Dᵢ = { 1 ; imag(λᵢ)>0
    #     { 0 ; imag(λᵢ)<0
    D = @. (sign(imag(λ))+1.0)/2
    Corr = V*Diagonal(D)*V⁻¹

    #Calculate the expected values for the currents
    aₖaₖ= diag(Corr)[1:L₁+L₂]
    aₖc = Corr[1:L₁+L₂, 1+2(L₁+L₂)]
    caₖ = Corr[1+2(L₁+L₂), 1:L₁+L₂]

    A = @. FermiDirac(ϵ, Vₗ, Tₗ)-real(aₖaₖ)
    B = @. real(aₖc)+real(caₖ)

    Jₚ = sum(@. γ*A)
    Jₕ = sum(@. γ*ϵ*A - √(Γ/(8π))*(γ^1.5)*B)

    return Jₚ, Jₕ
end

let
    #Machine parameters
    W  = 1.
    W° = W/2
    Γ₀ = 0.1W/8
    Tₗ = 1.1W/8
    Tᵣ = W/8
    L  = 100
    L₁ = Int(0.8*L)
    L₂ = Int(0.2*L)

    # Data frame to save the data
    result_df = DataFrame(ΔV = Float64[], μϵ_W = Float64[], η_ηc = Float64[], P_w2 = Float64[])

    #=
    # This is an equivalent way of doing it. We found that as long L is big enougth, they are the same. This one is a bit less stable.
    ϵ = 1/8
    # To do the 2D heat maps in Fig 8 of the Guide paper
    for ΔV in range(-0.1*W, 0.1*W, length = 50)
        println("ΔV = ", ΔV)
        for μ in range((-1.0 + ϵ)*W, (1.0 + ϵ)*W, length = 50)
            # Vᵣ + Vₗ = 2μ
            # Vᵣ - Vₗ = ΔV
            Vᵣ = (2μ + ΔV)/2
            Vₗ = (2μ - ΔV)/2
            Jₚ, Jₕ = RunMachine(ϵ, W, W°, Γ₀, Vₗ, Vᵣ, Tₗ, Tᵣ, L₁, L₂)

            dQₗ = Jₕ - Vₗ*Jₚ
            dQᵣ = Jₕ - Vᵣ*Jₚ
            P   = dQₗ - dQᵣ 
            η   = P/dQₗ
            η_c = 1.0 - Tᵣ/Tₗ

            # Save in the data frame:
            if 0 <= η && η <= η_c 
                push!(result_df, (ΔV, (μ - ϵ)/W, η/η_c, P/W^2))
            else
                push!(result_df, (ΔV, (μ - ϵ)/W, 0, 0))
            end

        end
    end
    =#

    # To do the 2D heat maps in Fig 8 of the Guide paper
    for ΔV in range(-0.1*W, 0.1*W, length = 80)
        println("ΔV = ", ΔV)
        for ϵ in range(-1.0*W, 1.0*W, length = 80)
            Vₗ = ΔV/2
            Vᵣ = -ΔV/2
            Jₚ, Jₕ = RunMachine(ϵ, W, W°, Γ₀, Vₗ, Vᵣ, Tₗ, Tᵣ, L₁, L₂)

            μ   = Vᵣ + Vₗ
            dQₗ = Jₕ - Vₗ*Jₚ
            dQᵣ = Jₕ - Vᵣ*Jₚ
            P   = dQₗ - dQᵣ 
            η   = 1.0 - dQᵣ/dQₗ
            η_c = 1.0 - Tᵣ/Tₗ

            # Save in the data frame:
            if 0 <= η && η <= η_c 
                push!(result_df, (ΔV, (μ - ϵ)/W, η/η_c, P/W^2))
            else
                push!(result_df, (ΔV, (μ - ϵ)/W, 0, 0))
            end

        end
    end

    # Save the data frame to a CSV file
    CSV.write("results/heatmap_data.csv", result_df)
    
    return
end
