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
    InnerRange = 1:L₁÷2
    Δϵ = 2W°/(L₁-1) 
    ϵ₁ = @. Δϵ*(InnerRange-0.5)
    γ₁ = @. Δϵ+0InnerRange

    # Logarithmic energies
    OuterRange = 1:L₂÷2
    Φ = ExponentialRate(L₂÷2, (W-W°)/Δϵ)
    ϵ₂ = @. W° + ((W-W°)/(1-Φ^(L₂÷2)))*(1.0-Φ^OuterRange)
    γ₂ = @. Δϵ*Φ.^OuterRange

    # Organice all
    ϵ = vcat(reverse(-ϵ₂), reverse(-ϵ₁), ϵ₁, ϵ₂)
    γ = vcat(reverse(γ₂), reverse(γ₁), γ₁, γ₂)

    return ϵ, γ
end
