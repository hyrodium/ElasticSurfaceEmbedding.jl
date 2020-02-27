using JSON

f=open("Julia/config.json", "r")
dict=JSON.parse(f)  # parse and transform data
const NIP=dict["Default Number of Integration Points"]
const 𝝂=dict["Poisson's Ratio"]
const OUT_DIR=dict["Output Directory"]
const IWhU = dict["Slack"]["Incoming Webhook Url"]
const OAAT = dict["Slack"]["OAuth Access Token"]
const ChID = dict["Slack"]["Channel ID"]
close(f)

const d=2 # Dimension
const Y=1.0 # Young Modulus
const 𝝀=𝝂*Y/((1+𝝂)*(1-(d-1)*𝝂)) # Lamé constant
const 𝝁=1/2(1+𝝂) # Lamé constant
