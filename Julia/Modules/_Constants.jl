f = open("Julia/config.json", "r")
config_dict = JSON.parse(f)  # parse and transform data
const NIP = config_dict["Default Number of Integration Points"]::Int
const 𝝂 = config_dict["Poisson's Ratio"]::Float64
const OUT_DIR = expanduser(config_dict["Output Directory"]::String)
const IWhU = config_dict["Slack"]["Incoming Webhook Url"]::String
const OAAT = config_dict["Slack"]["OAuth Access Token"]::String
const ChID = config_dict["Slack"]["Channel ID"]::String
close(f)

const d = 2 # Dimension
const Y = 1.0 # Young Modulus
const 𝝀 = 𝝂 * Y / ((1 + 𝝂) * (1 - (d - 1) * 𝝂)) # Lamé constant
const 𝝁 = 1 / 2(1 + 𝝂) # Lamé constant

const distributed = (@isdefined Distributed) && (Distributed isa Module)
const ESE_VERSION = v"0.0.1"
