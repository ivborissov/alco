using HetaSimulator, Plots

p = load_platform(".")

scenarios = read_scenarios("julia/scenarios.csv")
add_scenarios!(p, scenarios)



sim(p) |> plot
