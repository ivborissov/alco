using HetaSimulator, Plots

p = load_platform(".")

scenarios = read_scenarios("julia/scenarios.csv")
add_scenarios!(p, scenarios)

data = read_measurements("julia/data.csv")
add_measurements!(p, data)

sim(p) |> plot
