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
