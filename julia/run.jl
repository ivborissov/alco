using HetaSimulator, Plots

p = load_platform(".", rm_out=false)

scenarios = read_scenarios("data-mumenthaler-2000/scenarios.csv")
add_scenarios!(p, scenarios)

data = read_measurements("data-mumenthaler-2000/data.csv")
data_scn = Dict()
data_scn[:scn1] = filter(:scenario => ==(:scn1), data)
data_scn[:scn2] = filter(:scenario => ==(:scn2), data)
data_scn[:scn3] = filter(:scenario => ==(:scn3), data)

loss_add = Dict()
loss_add[:scn1] = 2*sum(log.(data_scn[:scn1][:,"prob.sigma"])) + length(data_scn[:scn1][:,"prob.sigma"])*log(2π)
loss_add[:scn2] = 2*sum(log.(data_scn[:scn2][:,"prob.sigma"])) + length(data_scn[:scn1][:,"prob.sigma"])*log(2π)
loss_add[:scn3] = 2*sum(log.(data_scn[:scn3][:,"prob.sigma"])) + length(data_scn[:scn1][:,"prob.sigma"])*log(2π)

add_measurements!(p, data)

# initial plot
sim(p) |> plot
sim(p, parameters_upd = [:Vmax=>0.3, :ke=>0., :k_a=>5.]) |> plot

# fitting 1 best -63.45
params_df_1 = read_parameters("./julia/parameters-1.csv")
res_fit_1 = fit(p, params_df_1; scenarios=[:scn1], ftol_abs=1e-6, ftol_rel=0.)
res_fit_1 = fit(p, optim(res_fit_1); scenarios=[:scn1], ftol_abs=1e-6, ftol_rel=0.)
fig = sim(p, scenarios=[:scn1], parameters_upd=optim(res_fit_1)) |> plot
# savefig(fig, "diagnostics/scn1_best.png")


# fitting 2 best -63.64
params_df_2 = read_parameters("./julia/parameters-2.csv")
res_fit_2 = fit(p, params_df_2; scenarios=[:scn2], ftol_abs=1e-6, ftol_rel=0.)
res_fit_2 = fit(p, optim(res_fit_2); scenarios=[:scn2], ftol_abs=1e-6, ftol_rel=0.)
fig = sim(p, scenarios=[:scn2], parameters_upd=optim(res_fit_2)) |> plot
# savefig(fig, "diagnostics/scn2_best.png")


# fitting 3 best -63.92
params_df_3 = read_parameters("./julia/parameters-3.csv")
res_fit_3 = fit(p, params_df_3; scenarios=[:scn3], ftol_abs=1e-6, ftol_rel=0.)
res_fit_3 = fit(p, optim(res_fit_3); scenarios=[:scn3], ftol_abs=1e-6, ftol_rel=0.)
fig = sim(p, scenarios=[:scn3], parameters_upd=optim(res_fit_3)) |> plot
# savefig(fig, "diagnostics/scn3_best.png")

############################## Identification ##########################

using LikelihoodProfiler, CSV

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
  sim_res = sim(p.scenarios[scen]; parameters_upd=params, reltol=1e-6, abstol=1e-8)
  #sim_res = last(sim_vec[1])
  return loss(sim_res, sim_res.scenario.measurements) - loss_add[scen]
end

function loss_func(pvals::Vector{N}, scen) where N <: Number
 # @show pvals
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
  theta_bounds = fill((1e-10,1e10), length(p_names[:scn1])),
  scan_bounds=((last.(p_optim_1)[i])/1e4, (last.(p_optim_1)[i])*1e4),
  scan_tol=1e-5,
  #scale =  fill(:log, length(p_names[:scn1])),
  loss_crit = loss_scn1(p_optim_1) + chi_level
  ) for i in eachindex(p_names[:scn1])]

plb_1 = [iv.result[1].value for iv in p_ident_1]
pub_1 = [iv.result[2].value for iv in p_ident_1]

pliter_1 = [iv.result[1].counter for iv in p_ident_1]
puiter_1 = [iv.result[2].counter for iv in p_ident_1]

df = DataFrame(params = p_names[:scn1], optim = last.(p_optim_1), lower = plb_1, upper = pub_1, liter = pliter_1, uiter = puiter_1)
df.lower = replace(df.lower, nothing => missing)
df.upper = replace(df.upper, nothing => missing)
CSV.write("./julia/scn1_intervals.csv", df)

p_ident_2 = [get_interval(
  last.(p_optim_2),
  i,
  loss_scn2,
  :CICO_ONE_PASS,
  #theta_bounds = fill((-10.,10.), length(p_names[:scn2])),
  #scale = fill(:log, length(p_names[:scn2])),
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

liter_1 = [iv.result[1].counter for iv in BrAC_ident_1]
uiter_1 = [iv.result[2].counter for iv in BrAC_ident_1]

df = DataFrame(times = saveat_1, lower = lb_1, upper = ub_1, liter = liter_1, uiter = uiter_1)
df.lower = replace(df.lower, nothing => missing)
df.upper = replace(df.upper, nothing => missing)
CSV.write("./julia/scn1_conf_band.csv", df)

BrAC_1 = sim_scn1.(saveat_1, :BrAC)
plot(sim_scn1, show_measurements=false, vars=[:BrAC])
scatter!(data_scn[:scn1].t, data_scn[:scn1].measurement, yerror=data_scn[:scn1][!,"prob.sigma"], label = "Measurements")
plot!(saveat_1, lb_1, fillrange =  ub_1, fillalpha = 0.35, c = 1, label = "Confidence band")
savefig("./julia/conf_band.png")

lb_2 = [iv.result[1].value for iv in BrAC_ident_2]
ub_2 = [iv.result[2].value for iv in BrAC_ident_2]
BrAC_2 = sim_scn2.(saveat_2, :BrAC)
plot(sim_scn2, show_measurements=false, ribbon = (BrAC_2-lb_1,ub_2-BrAC_2), fc=:orange, fa=0.7)

lb_3 = [iv.result[1].value for iv in BrAC_ident_3]
ub_3 = [iv.result[2].value for iv in BrAC_ident_3]
BrAC_3 = sim_scn3.(saveat_3, :BrAC)
plot(sim_scn3, show_measurements=false, ribbon = (BrAC_3-lb_3,ub_3-BrAC_3), fc=:orange, fa=0.7)

### Validation band

t_scn1 = data_scn[:scn1].t
sigma_scn1 = data_scn[:scn1][!,"prob.sigma"]
sim_scn1 = sim(p.scenarios[:scn1]; parameters_upd=p_optim_1, reltol=1e-6, abstol=1e-8)
BrAC_scn1 = sim_scn1.(t_scn1, :BrAC)

function valid_obj1(params, i)
  d1 = last(params)
  _params = [pn => pv for (pn,pv) in zip(p_names[:scn1],params[1:end-1])]

  sim_res = sim(p.scenarios[:scn1]; parameters_upd=_params, reltol=1e-6, abstol=1e-8)
  d1_sim = sim_res(t_scn1[i])[:BrAC]
  return loss(sim_res, sim_res.scenario.measurements) + (d1 - d1_sim)^2/(sigma_scn1[i])^2 - loss_add[:scn1]
end

valid_1 = []
for i in eachindex(t_scn1) 
println(" Calculating CI for $(t_scn1[i]) timepoint")
push!(valid_1, 
  get_interval(
    [last.(p_optim_1); BrAC_scn1[i]],
    5,
    p->valid_obj1(p,i),
    :CICO_ONE_PASS,
    theta_bounds = fill((1e-8,1e8), length(p_names[:scn1])+1),
    scan_bounds=(1e-7,1e7), 
    scan_tol=1e-5,
    scale =  fill(:log, length(p_names[:scn1])+1),
    loss_crit = loss_scn1(p_optim_1) + chi_level)
  )
end

lb_1 = [iv.result[1].value for iv in valid_1]
ub_1 = [iv.result[2].value for iv in valid_1]

liter_1 = [iv.result[1].counter for iv in valid_1]
uiter_1 = [iv.result[2].counter for iv in valid_1]

df = DataFrame(times = t_scn1, lower = lb_1, upper = ub_1, liter = liter_1, uiter = uiter_1)
df.lower = replace(df.lower, nothing => missing)
df.upper = replace(df.upper, nothing => missing)
CSV.write("./julia/scn1_valid_bans.csv", df)

lb_1[1] = 0.0
lb_1[end-4:end] .= 0.0
plot(sim_scn1, show_measurements=false, vars=[:BrAC])
scatter!(data_scn[:scn1].t, data_scn[:scn1].measurement, yerror=data_scn[:scn1][!,"prob.sigma"], label = "Measurements")
plot!(t_scn1, lb_1, fillrange =  ub_1, fillalpha = 0.35, c = 1, label = "Validation band")
savefig("./julia/valid_band.png")