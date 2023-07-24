## Reference State

# Parametric mapping of Reference state
𝒑₍₀₎(u¹, u²) = SVector(u¹, u², 0)
const surface = 𝒑₍₀₎

# Tangent vector
𝒑₁₍₀₎(u¹, u²) = ForwardDiff.derivative(u¹ -> 𝒑₍₀₎(u¹, u²), u¹)
𝒑₂₍₀₎(u¹, u²) = ForwardDiff.derivative(u² -> 𝒑₍₀₎(u¹, u²), u²)
𝒑₁₁₍₀₎(u¹, u²) = ForwardDiff.derivative(u¹ -> 𝒑₁₍₀₎(u¹, u²), u¹)
𝒑₁₂₍₀₎(u¹, u²) = ForwardDiff.derivative(u² -> 𝒑₁₍₀₎(u¹, u²), u²)
𝒑₂₁₍₀₎(u¹, u²) = ForwardDiff.derivative(u¹ -> 𝒑₂₍₀₎(u¹, u²), u¹)
𝒑₂₂₍₀₎(u¹, u²) = ForwardDiff.derivative(u² -> 𝒑₂₍₀₎(u¹, u²), u²)

# Normal vector
𝒆₍₀₎(u¹, u²) = normalize(cross(𝒑₁₍₀₎(u¹, u²), 𝒑₂₍₀₎(u¹, u²)))

# Riemannian metrix
g₍₀₎₁₁(u¹, u²) = dot(𝒑₁₍₀₎(u¹, u²), 𝒑₁₍₀₎(u¹, u²))
g₍₀₎₁₂(u¹, u²) = dot(𝒑₁₍₀₎(u¹, u²), 𝒑₂₍₀₎(u¹, u²))
g₍₀₎₂₁(u¹, u²) = dot(𝒑₂₍₀₎(u¹, u²), 𝒑₁₍₀₎(u¹, u²))
g₍₀₎₂₂(u¹, u²) = dot(𝒑₂₍₀₎(u¹, u²), 𝒑₂₍₀₎(u¹, u²))
g₍₀₎(u¹, u²) = @SMatrix [g₍₀₎₁₁(u¹, u²) g₍₀₎₁₂(u¹, u²); g₍₀₎₂₁(u¹, u²) g₍₀₎₂₂(u¹, u²)]
h₍₀₎(u¹, u²) = @SMatrix [
    (𝒆₍₀₎(u¹, u²)'*𝒑₁₁₍₀₎(u¹, u²)) (𝒆₍₀₎(u¹, u²)'*𝒑₁₂₍₀₎(u¹, u²))
    (𝒆₍₀₎(u¹, u²)'*𝒑₂₁₍₀₎(u¹, u²)) (𝒆₍₀₎(u¹, u²)'*𝒑₂₂₍₀₎(u¹, u²))
]

# Gaussian curvature
K₍₀₎(u¹, u²) = det(h₍₀₎(u¹, u²)) / det(g₍₀₎(u¹, u²))

# Volume form
𝝊₍₀₎(u¹, u²) = norm(cross(𝒑₁₍₀₎(u¹, u²), 𝒑₂₍₀₎(u¹, u²)))
g⁻₍₀₎(u¹, u²) = inv(g₍₀₎(u¹, u²)) # 第1基本量の逆
g₁₍₀₎(u¹, u²) = ForwardDiff.derivative(u¹ -> g₍₀₎(u¹, u²), u¹)
g₂₍₀₎(u¹, u²) = ForwardDiff.derivative(u² -> g₍₀₎(u¹, u²), u²)

# Christoffel symbol
𝛤₍₀₎²₁₁(u¹, u²) =
    (g⁻₍₀₎(u¹, u²)[2, 1] * g₁₍₀₎(u¹, u²)[1, 1] + g⁻₍₀₎(u¹, u²)[2, 2] * (2g₁₍₀₎(u¹, u²)[2, 1] - g₂₍₀₎(u¹, u²)[1, 1])) / 2
e⁽⁰⁾₁(u¹, u²) = normalize(𝒑₁₍₀₎(u¹, u²))
e⁽⁰⁾₂(u¹, u²) = normalize(𝒑₂₍₀₎(u¹, u²) - (g₍₀₎₁₂(u¹, u²) / g₍₀₎₁₁(u¹, u²)) * 𝒑₁₍₀₎(u¹, u²))

c(D₂::ClosedInterval) = sum(extrema(D₂)) / 2 # Coordinate on the center curve
s₍₀₎(t, D₂::ClosedInterval) = sqrt(g₍₀₎₁₁(t, c(D₂)))
ṡ₍₀₎(t, D₂::ClosedInterval) = (1 / 2) * (g₁₍₀₎(t, c(D₂)))[1, 1] / sqrt(g₍₀₎₁₁(t, c(D₂)))
𝜅₍₀₎(t, D₂::ClosedInterval) = 𝛤₍₀₎²₁₁(t, c(D₂)) * 𝝊₍₀₎(t, c(D₂)) / s₍₀₎(t, D₂)^3 # Geodesic curvature
K₍₀₎(t, D₂::ClosedInterval) = K₍₀₎(t, c(D₂)) # Gaussian curvature
B̃(t, D₂::ClosedInterval) = dot(e⁽⁰⁾₂(t, c(D₂)), 𝒑₂₍₀₎(t, c(D₂))) * width(D₂) / 2 # Breadth of the piece of surface
g₍₀₎₁₁(u¹, D₂::ClosedInterval) = g₍₀₎₁₁(u¹, c(D₂))
g₍₀₎₁₂(u¹, D₂::ClosedInterval) = g₍₀₎₁₂(u¹, c(D₂))
g₍₀₎₂₁(u¹, D₂::ClosedInterval) = g₍₀₎₂₁(u¹, c(D₂))
g₍₀₎₂₂(u¹, D₂::ClosedInterval) = g₍₀₎₂₂(u¹, c(D₂))
𝝊₍₀₎(u¹, D₂::ClosedInterval) = 𝝊₍₀₎(u¹, c(D₂))

# Current State
𝒑₍ₜ₎(M, u¹, u²) = M(u¹, u²)
# This can be faster with BSplineDerivativeSpace, but we don't need speed here.
𝒑₁₍ₜ₎(M, u¹, u²) = ForwardDiff.derivative(u¹ -> 𝒑₍ₜ₎(M, u¹, u²), u¹)
𝒑₂₍ₜ₎(M, u¹, u²) = ForwardDiff.derivative(u² -> 𝒑₍ₜ₎(M, u¹, u²), u²)

g₍ₜ₎₁₁(M, u¹, u²) = dot(𝒑₁₍ₜ₎(M, u¹, u²), 𝒑₁₍ₜ₎(M, u¹, u²)) # 第1基本量
g₍ₜ₎₁₂(M, u¹, u²) = dot(𝒑₁₍ₜ₎(M, u¹, u²), 𝒑₂₍ₜ₎(M, u¹, u²)) # 第1基本量
g₍ₜ₎₂₁(M, u¹, u²) = dot(𝒑₂₍ₜ₎(M, u¹, u²), 𝒑₁₍ₜ₎(M, u¹, u²)) # 第1基本量
g₍ₜ₎₂₂(M, u¹, u²) = dot(𝒑₂₍ₜ₎(M, u¹, u²), 𝒑₂₍ₜ₎(M, u¹, u²)) # 第1基本量
g₍ₜ₎(M, u¹, u²) = @SMatrix [g₍ₜ₎₁₁(M, u¹, u²) g₍ₜ₎₁₂(M, u¹, u²); g₍ₜ₎₂₁(M, u¹, u²) g₍ₜ₎₂₂(M, u¹, u²)]
