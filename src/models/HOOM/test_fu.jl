using Formatting

α_AB = 1.0 /  2.0
β_AB = 5.0 / 12.0
@inline function ABIII(a, aa, aaa)
    return ( 1 + α_AB + β_AB ) * a - (α_AB + 2*β_AB) * aa + β_AB * aaa
end
 


τ = 1e-1
f = 1e-3
x = [ 1.0, 0 ]

α = 1.5


Δt = τ / f

all_steps = floor(10 / τ)

A = f * reshape([ 0, -1.0, 1.0, 0.0], 2,2)

G = zeros(Float64, 2, 3)

println(format("[00] x = [ {:.2f}, {:.2f} ] ", x...))
for t = 1:all_steps
    
    G[:, 1] .= A * x
    
    #G_aux = α * G[:, 1] + ( 1 - α ) * G[:, 2]
    G_aux = ABIII(G[:, 1], G[:, 2], G[:, 3])

    x .+= Δt * G_aux

    G[:, 2] = G[:, 1]
    G[:, 3] = G[:, 2]
    
    println(format("[{:02d}] x = [ {:.2f}, {:.2f} ] ; length = {:.2f}", t, x..., sum(x.^2)^.5))

end





