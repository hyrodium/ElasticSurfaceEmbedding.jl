const d = 2 # Dimension of surface

const 𝝂 = 0.25 # Poisson's Ratio
const Y = 1.0 # Young Modulus
const 𝝀 = 𝝂 * Y / ((1 + 𝝂) * (1 - (d - 1) * 𝝂)) # Lamé constant
const 𝝁 = 1 / 2(1 + 𝝂) # Lamé constant
