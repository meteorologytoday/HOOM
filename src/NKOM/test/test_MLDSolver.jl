include("../NKOM.jl")
using .NKOM


h = range(0, 100, length=1000) |> collect

Δ    = copy(h)
∂Δ∂h = copy(h)


Bf = 330.0
J0 = -1000.0
fric_u = 0.1
f = 1e-4
m = 0.8
n = 0.2
γ = 1/23.0

a, λ = NKOM.calΔCoefficients(fric_u, f, m)

for i = 1:length(h)
    Δ[i], ∂Δ∂h[i] = NKOM.calΔ_and_∂Δ∂h(h[i], a, Bf, J0, λ, γ, n)
end

using PyCall

plt = pyimport("matplotlib.pyplot")


fig, ax = plt.subplots(2, 1, sharex=true)

ax[1].plot(h, Δ, label="Δ")
ax[2].plot(h, ∂Δ∂h, label="∂Δ∂h")

ax[1].legend()
ax[2].legend()

plt.show()
