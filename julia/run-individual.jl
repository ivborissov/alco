using HetaSimulator, Plots

p = load_platform(".")

scenarios = read_scenarios("data-mumenthaler-2000/scenarios.csv")
add_scenarios!(p, scenarios)

data = read_measurements("data-mumenthaler-2000/individual.csv")
add_measurements!(p, data)

# initial plot
sim(p) |> plot
sim(p, parameters_upd = [:Vmax=>0.3, :ke=>0., :k_a=>5.]) |> plot

# fitting 1 best -73.48/-80.87, sigma_add = 0.035/ sigma_prop=0.081
params_df_1 = read_parameters("./julia/parameters-1-ind.csv")
res_fit_1 = fit(p, params_df_1; scenarios=[:scn1], ftol_abs=1e-6, ftol_rel=0.)
res_fit_1 = fit(p, optim(res_fit_1); scenarios=[:scn1], ftol_abs=1e-6, ftol_rel=0.)
fig = sim(p, scenarios=[:scn1], parameters_upd=optim(res_fit_1)) |> plot
# savefig(fig, "diagnostics/scn1_best.png")
