using HetaSimulator, Plots

p = load_platform(".")

scenarios = read_scenarios("julia/scenarios.csv")
add_scenarios!(p, scenarios)

data = read_measurements("julia/data.csv")
add_measurements!(p, data)

# initial plot
sim(p) |> plot
sim(p, parameters_upd = [:Vmax=>0.3, :ke=>0., :k_a=>5.]) |> plot

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


# fitting 3 best -126.25 / -95.83
params_df_3 = read_parameters("./julia/parameters-3.csv")
res_fit_3 = fit(p, params_df_3; scenarios=[:scn3])
res_fit_3 = fit(p, optim(res_fit_3); scenarios=[:scn3])
fig = sim(p, scenarios=[:scn3], parameters_upd=optim(res_fit_3)) |> plot
# savefig(fig, "diagnostics/scn3_best.png")

############################## Identification ##########################

using LikelihoodProfiler

chi_level = 3.84

p_optim_1 = optim(res_fit_1)
p_optim_2 = optim(res_fit_2)
p_optim_3 = optim(res_fit_3)

sim_scn1 = sim(p.scenarios[:scn1], parameters_upd=p_optim_1)
sim_scn2 = sim(p.scenarios[:scn2], parameters_upd=p_optim_2)
sim_scn3 = sim(p.scenarios[:scn3], parameters_upd=p_optim_3)

p_names = Dict(
  :scn1 => first.(p_optim_1),
  :scn2 => first.(p_optim_2),
  :scn3 => first.(p_optim_3),
)

function loss_func(params::Vector{P}, scen) where P<: Pair
  sim_vec = sim(p; parameters_upd=params, scenarios=[scen])
  sim_res = last(sim_vec[1])
  return loss(sim_res, sim_res.scenario.measurements)
end

function loss_func(pvals::Vector{N}, scen) where N <: Number
  @assert length(p_names[scen]) == length(pvals) "Number of params values doesn't match params names"
  params = [pn => pv for (pn,pv) in zip(p_names[scen],pvals)]
  loss_func(params,scen)
end

loss_scn1(params) = loss_func(params, :scn1)
loss_scn2(params) = loss_func(params, :scn2)
loss_scn3(params) = loss_func(params, :scn3)

function scan_func(params::Vector{P}, timepoint, scen) where P<: Pair
  sim_vec = sim(p; parameters_upd=params, scenarios=[scen])
  return last(sim_vec[1])(timepoint)[:BrAC]
end

function scan_func(pvals::Vector{N}, timepoint, scen) where N <: Number
  @assert length(p_names[scen]) == length(pvals) "Number of params values doesn't match params names"
  params = [pn => pv for (pn,pv) in zip(p_names[scen],pvals)]
  scan_func(params, timepoint, scen)
end

scan_scn1(params, timepoint) = scan_func(params, timepoint, :scn1)
scan_scn2(params, timepoint) = scan_func(params, timepoint, :scn2)
scan_scn3(params, timepoint) = scan_func(params, timepoint, :scn3)

saveat_1 = saveat(p.scenarios[:scn1])
saveat_2 = saveat(p.scenarios[:scn2])
saveat_3 = saveat(p.scenarios[:scn3])

p_ident_1 = [get_interval(
  last.(p_optim_1),
  i,
  loss_scn1,
  :CICO_ONE_PASS,
  scale = fill(:log, length(p_names[:scn1])),
  loss_crit = loss_scn1(p_optim_1) + chi_level
  ) for i in eachindex(p_names[:scn1])]

p_ident_2 = [get_interval(
  last.(p_optim_2),
  i,
  loss_scn2,
  :CICO_ONE_PASS,
  scale = fill(:log, length(p_names[:scn2])),
  loss_crit = loss_scn2(p_optim_2) + chi_level
) for i in eachindex(p_names[:scn2])]

p_ident_3 = [get_interval(
  last.(p_optim_3),
  i,
  loss_scn3,
  :CICO_ONE_PASS,
  scale = fill(:log, length(p_names[:scn3])),
  loss_crit = loss_scn3(p_optim_3) + chi_level
) for i in eachindex(p_names[:scn3])]

BrAC_ident_1 = [get_interval(
  last.(p_optim_1),
  params->scan_scn1(params,t),
  loss_scn1,
  :CICO_ONE_PASS,
  scale = fill(:log, length(p_names[:scn1])),
  loss_crit = loss_scn1(p_optim_1) + chi_level
) for t in saveat_1]

BrAC_ident_2 = [get_interval(
  last.(p_optim_2),
  params->scan_scn2(params,t),
  loss_scn2,
  :CICO_ONE_PASS,
  scale = fill(:log, length(p_names[:scn2])),
  loss_crit = loss_scn2(p_optim_2) + chi_level
) for t in saveat_2]

BrAC_ident_3 = [get_interval(
  last.(p_optim_3),
  params->scan_scn3(params,t),
  loss_scn3,
  :CICO_ONE_PASS,
  scale = fill(:log, length(p_names[:scn3])),
  loss_crit = loss_scn3(p_optim_3) + chi_level
) for t in saveat_3]

lb_1 = [iv.result[1].value for iv in BrAC_ident_1]
ub_1 = [iv.result[2].value for iv in BrAC_ident_1]
BrAC_1 = sim_scn1.(saveat_1, :BrAC)
plot(sim_scn1, show_measurements=false, ribbon = (BrAC_1-lb_1,ub_1-BrAC_1), fc=:orange, fa=0.7)
#savefig("BrAC_ident_scn1.png")

lb_2 = [iv.result[1].value for iv in BrAC_ident_2]
ub_2 = [iv.result[2].value for iv in BrAC_ident_2]
BrAC_2 = sim_scn2.(saveat_2, :BrAC)
plot(sim_scn2, show_measurements=false, ribbon = (BrAC_2-lb_1,ub_2-BrAC_2), fc=:orange, fa=0.7)

lb_3 = [iv.result[1].value for iv in BrAC_ident_3]
ub_3 = [iv.result[2].value for iv in BrAC_ident_3]
BrAC_3 = sim_scn3.(saveat_3, :BrAC)
plot(sim_scn3, show_measurements=false, ribbon = (BrAC_3-lb_3,ub_3-BrAC_3), fc=:orange, fa=0.7)