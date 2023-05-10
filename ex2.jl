using ITensors
N = 3

ITensors.op(::OpName"σx", ::SiteType"S=1/2", s::Index) =
  2*op("Sx", s)

function ITensors.op(::OpName"σz", ::SiteType"S=1/2", s::Index)
    mat = [1.0 0.0
           0.0 -1.0]
    return itensor(mat, s', s)
end
#2*op("Sz", s)

# Make the operator list.
os = [("σx", n) for n in 1:N]
append!(os, [("σz", n) for n in 1:N])

@show os

s = siteinds("S=1/2", N)
gates = ops(os, s)

# Starting state |↑↑↑⟩
ψ0 = MPS(s, "↑")

# Apply the gates.
ψ = apply(gates, ψ0; cutoff = 1e-15)

# Test against exact (full) wavefunction
prodψ = apply(gates, prod(ψ0))
@show prod(ψ) ≈ prodψ

# The result is:
# σz₃ σz₂ σz₁ σx₃ σx₂ σx₁ |↑↑↑⟩ = -|↓↓↓⟩
@show inner(ψ, MPS(s, "↓")) == -1

# 2-site gate
function ITensors.op(::OpName"CX", ::SiteType"S=1/2", s1::Index, s2::Index)
  mat = [1 0 0 0
         0 1 0 0
         0 0 0 1
         0 0 1 0]
  return itensor(mat, s2', s1', s2, s1)
end

os = [("CX", 1, 3), ("σz", 3)]

@show os
@show op("CX", s[1], s[2]) - op("σz * σz", s[1])*op("σz * σz", s[2])

# Start with the state |↓↑↑⟩
ψ0 = MPS(s, n -> n == 1 ? "↓" : "↑")

# The result is:
# σz₃ CX₁₃ |↓↑↑⟩ = -|↓↑↓⟩
ψ = apply(ops(os, s), ψ0; cutoff = 1e-15)
@show inner(ψ, MPS(s, n -> n == 1 || n == 3 ? "↓" : "↑")) == -1
