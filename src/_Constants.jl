const d = 2 # Dimension of surface
const NIP = 10 # Default number of integration points for Gaussan quadrature

const 𝝂 = 0.25 # Poisson's Ratio
const Y = 1.0 # Young Modulus
const 𝝀 = 𝝂 * Y / ((1 + 𝝂) * (1 - (d - 1) * 𝝂)) # Lamé constant
const 𝝁 = 1 / 2(1 + 𝝂) # Lamé constant

const ESE_VERSION = v"0.0.1"

# Default output directory
OUT_DIR = joinpath(homedir(),"ElasticSurfaceEmbedding-Result")

"""
    config_dir(dir)

Set the output directory.
The default is `~/ElasticSurfaceEmbedding-Result`.
"""
function config_dir(dir)
    _dir = expanduser(dir)
    mkpath(_dir)
    global OUT_DIR = _dir
end

struct SlackConfig
    channel::String
    token::String
end

SLACK = SlackConfig("","")

"""
config_slack(;channel, token)
"""
function config_slack(;channel, token)
    global SLACK = SlackConfig(channel, token)
end
