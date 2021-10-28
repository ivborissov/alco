using HetaSimulator, Plots

p = load_platform(".")

scenarios = read_scenarios("julia/scenarios.csv")
add_scenarios!(p, scenarios)

data = read_measurements("julia/data.csv")
add_measurements!(p, data)

# initial plot
sim(p) |> plot

# fitting 1
params_df_1 = read_parameters("./julia/parameters-1.csv")
res_fit_1 = fit(p, params_df_1; scenarios=[:scn1])
res_fit_1 = fit(p, optim(res_fit_1); scenarios=[:scn1])
sim(p, parameters_upd=optim(res_fit_1)) |> plot
# best -96.12

# fitting 2
params_df_2 = read_parameters("./julia/parameters-2.csv")
res_fit_2 = fit(p, params_df_2; scenarios=[:scn2])
res_fit_2 = fit(p, optim(res_fit_2); scenarios=[:scn2])
sim(p, parameters_upd=optim(res_fit_2)) |> plot
# best -102.04

# fitting 3
params_df_3 = read_parameters("./julia/parameters-3.csv")
res_fit_3 = fit(p, params_df_3; scenarios=[:scn3])
res_fit_3 = fit(p, optim(res_fit_3); scenarios=[:scn3])
sim(p, parameters_upd=optim(res_fit_3)) |> plot
# best -126.25
