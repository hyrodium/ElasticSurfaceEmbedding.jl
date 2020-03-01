using JSON

f=open("Julia/config.json", "r")
dict=JSON.parse(f)  # parse and transform data
const NIP=dict["Default Number of Integration Points"] :: Int
const 𝝂=dict["Poisson's Ratio"] :: Float64
const OUT_DIR=expanduser(dict["Output Directory"] :: String)
const IWhU = dict["Slack"]["Incoming Webhook Url"] :: String
const OAAT = dict["Slack"]["OAuth Access Token"] :: String
const ChID = dict["Slack"]["Channel ID"] :: String
close(f)

const d=2 # Dimension
const Y=1.0 # Young Modulus
const 𝝀=𝝂*Y/((1+𝝂)*(1-(d-1)*𝝂)) # Lamé constant
const 𝝁=1/2(1+𝝂) # Lamé constant
