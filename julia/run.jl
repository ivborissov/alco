using HetaSimulator, Plots

p = load_platform(".")

scenarios = read_scenarios("julia/scenarios.csv")
add_scenarios!(p, scenarios)

data = read_measurements("julia/data.csv")
add_measurements!(p, data)

# initial plot
sim(p) |> plot

# fitting 1 best -96.12
params_df_1 = read_parameters("./julia/parameters-1.csv")
res_fit_1 = fit(p, params_df_1; scenarios=[:scn1])
res_fit_1 = fit(p, optim(res_fit_1); scenarios=[:scn1])
fig = sim(p, scenarios=[:scn1], parameters_upd=optim(res_fit_1)) |> plot
# savefig(fig, "diagnostics/scn1_best.png")


# fitting 2 best -102.04
params_df_2 = read_parameters("./julia/parameters-2.csv")
res_fit_2 = fit(p, params_df_2; scenarios=[:scn2])
res_fit_2 = fit(p, optim(res_fit_2); scenarios=[:scn2])
fig = sim(p, scenarios=[:scn2], parameters_upd=optim(res_fit_2)) |> plot
# savefig(fig, "diagnostics/scn2_best.png")


# fitting 3 best -126.25
params_df_3 = read_parameters("./julia/parameters-3.csv")
res_fit_3 = fit(p, params_df_3; scenarios=[:scn3])
res_fit_3 = fit(p, optim(res_fit_3); scenarios=[:scn3])
fig = sim(p, scenarios=[:scn3], parameters_upd=optim(res_fit_3)) |> plot
# savefig(fig, "diagnostics/scn3_best.png")

