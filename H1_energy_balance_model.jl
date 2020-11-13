### A Pluto.jl notebook ###
# v0.12.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1e06178a-1fbf-11eb-32b3-61769a79b7c0
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			"Plots",
			"PlutoUI",
			"LaTeXStrings",
			"Distributions",
			"Random",
	])
	using LaTeXStrings
	using Plots
	using PlutoUI
	using Random, Distributions
	Random.seed!(123)
	
	md"##### Package dependencies"
end

# ╔═╡ 169727be-2433-11eb-07ae-ab7976b5be90
md"_homework 9, version 0_"

# ╔═╡ 21524c08-2433-11eb-0c55-47b1bdc9e459
md"""

# **Homework 9**: _Climate modeling I_
`18.S191`, fall 2020

This notebook contains _built-in, live answer checks_! In some exercises you will see a coloured box, which runs a test case on your code, and provides feedback based on the result. Simply edit the code, run it, and the check runs again.

_For MIT students:_ there will also be some additional (secret) test cases that will be run as part of the grading process, and we will look at your notebook and write comments.

Feel free to ask questions!
"""

# ╔═╡ 23335418-2433-11eb-05e4-2b35dc6cca0e
# edit the code below to set your name and kerberos ID (i.e. email without @mit.edu)

student = (name = "Jazzy Doe", kerberos_id = "jazz")

# you might need to wait until all other cells in this notebook have completed running. 
# scroll around the page to see what's up

# ╔═╡ 18be4f7c-2433-11eb-33cb-8d90ca6f124c
md"""

Submission by: **_$(student.name)_** ($(student.kerberos_id)@mit.edu)
"""

# ╔═╡ 253f4da0-2433-11eb-1e48-4906059607d3
md"_Let's create a package environment:_"

# ╔═╡ 87e68a4a-2433-11eb-3e9d-21675850ed71
html"""
<iframe width="100%" height="300" src="https://www.youtube.com/embed/Gi4ZZVS2GLA" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 1312525c-1fc0-11eb-2756-5bc3101d2260
md"""## **Exercise 1** - _policy goals under uncertainty_
A recent ground-breaking [review paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019RG000678) produced the most comprehensive and up-to-date estimate of the *climate feedback parameter*, which they find to be

$B \approx \mathcal{N}(-1.3, 0.4),$

i.e. our knowledge of the real value is normally distributed with a mean value $\overline{B} = -1.3$ W/m²/K and a standard deviation $\sigma = 0.4$ W/m²/K. These values are not very intuitive, so let us convert them into more policy-relevant numbers.

**Definition:** *Equilibrium climate sensitivity (ECS)* is defined as the amount of warming $\Delta T$ caused by a doubling of CO₂ (e.g. from the pre-industrial value 280 ppm to 560 ppm), at equilibrium.

At equilibrium, the energy balance model equation is:

$0 = \frac{S(1 - α)}{4} - (A - BT_{eq}) + a \ln\left( \frac{2\;\text{CO}₂_{\text{PI}}}{\text{CO}₂_{\text{PI}}} \right)$

From this, we subtract the preindustrial energy balance, which is given by:

$0 = \frac{S(1-α)}{4} - (A - BT_{0}),$

The result of this subtraction, after rearranging, is our definition of $\text{ECS}$:

$\text{ECS} \equiv T_{eq} - T_{0} = -\frac{a\ln(2)}{B}$
"""

# ╔═╡ 7f961bc0-1fc5-11eb-1f18-612aeff0d8df
md"""The plot below provides an example of an "abrupt 2xCO₂" experiment, a classic experimental treatment method in climate modelling which is used in practice to estimate ECS for a particular model (Note: in complicated climate models the values of the parameters $a$ and $B$ are not specified *a priori*, but *emerge* as outputs of the simulation).

The simulation begins at the preindustrial equilibrium, i.e. a temperature $T_{0} = 14$°C is in balance with the pre-industrial CO₂ concentration of 280 ppm until CO₂ is abruptly doubled from 280 ppm to 560 ppm. The climate responds by rapidly warming, and after a few hundred years approaches the equilibrium climate sensitivity value, by definition.
"""

# ╔═╡ 930d7154-1fbf-11eb-1c3a-b1970d291811
module Model

const S = 1368; # solar insolation [W/m^2]  (energy per unit time per unit area)
const α = 0.3; # albedo, or planetary reflectivity [unitless]
const B = -1.3; # climate feedback parameter [W/m^2/°C],
const T0 = 14.; # preindustrial temperature [°C]

absorbed_solar_radiation(; α=α, S=S) = S*(1 - α)/4; # [W/m^2]
outgoing_thermal_radiation(T; A=A, B=B) = A - B*T;

const A = S*(1. - α)/4 + B*T0; # [W/m^2].

greenhouse_effect(CO2; a=a, CO2_PI=CO2_PI) = a*log(CO2/CO2_PI);

const a = 5.0; # CO2 forcing coefficient [W/m^2]
const CO2_PI = 280.; # preindustrial CO2 concentration [parts per million; ppm];
CO2_const(t) = CO2_PI; # constant CO2 concentrations

const C = 51.; # atmosphere and upper-ocean heat capacity [J/m^2/°C]

function timestep!(ebm)
	append!(ebm.T, ebm.T[end] + ebm.Δt*tendency(ebm));
	append!(ebm.t, ebm.t[end] + ebm.Δt);
end;

tendency(ebm) = (1. /ebm.C) * (
	+ absorbed_solar_radiation(α=ebm.α, S=ebm.S)
	- outgoing_thermal_radiation(ebm.T[end], A=ebm.A, B=ebm.B)
	+ greenhouse_effect(ebm.CO2(ebm.t[end]), a=ebm.a, CO2_PI=ebm.CO2_PI)
);

begin
	mutable struct EBM
		T::Array{Float64, 1}
	
		t::Array{Float64, 1}
		Δt::Float64
	
		CO2::Function
	
		C::Float64
		a::Float64
		A::Float64
		B::Float64
		CO2_PI::Float64
	
		α::Float64
		S::Float64
	end;
	
	# Make constant parameters optional kwargs
	EBM(T::Array{Float64, 1}, t::Array{Float64, 1}, Δt::Float64, CO2::Function;
		C=C, a=a, A=A, B=B, CO2_PI=CO2_PI, α=α, S=S) = (
		EBM(T, t, Δt, CO2, C, a, A, B, CO2_PI, α, S)
	);
	
	# Construct from float inputs for convenience
	EBM(T0::Float64, t0::Float64, Δt::Float64, CO2::Function;
		C=C, a=a, A=A, B=B, CO2_PI=CO2_PI, α=α, S=S) = (
		EBM([T0], [t0], Δt, CO2;
			C=C, a=a, A=A, B=B, CO2_PI=CO2_PI, α=α, S=S);
	);
end;

begin
	function run!(ebm::EBM, end_year::Real)
		while ebm.t[end] < end_year
			timestep!(ebm)
		end
	end;
	
	run!(ebm) = run!(ebm, 200.) # run for 200 years by default
end




CO2_hist(t) = CO2_PI * (1 .+ fractional_increase(t));
fractional_increase(t) = ((t .- 1850.)/220).^3;

begin
	CO2_RCP26(t) = CO2_PI * (1 .+ fractional_increase(t) .* min.(1., exp.(-((t .-1850.).-170)/100))) ;
	RCP26 = EBM(T0, 1850., 1., CO2_RCP26)
	run!(RCP26, 2100.)
	
	CO2_RCP85(t) = CO2_PI * (1 .+ fractional_increase(t) .* max.(1., exp.(((t .-1850.).-170)/100)));
	RCP85 = EBM(T0, 1850., 1., CO2_RCP85)
	run!(RCP85, 2100.)
end

end

# ╔═╡ fa7e6f7e-2434-11eb-1e61-1b1858bb0988
md"""
``B = `` $(@bind B_slider Slider(-2.5:.001:0; show_value=true, default=-1.3))
"""

# ╔═╡ 16348b6a-1fc2-11eb-0b9c-65df528db2a1
md"""
##### Exercise 1.1 - _Develop understanding for feedbacks and climate sensitivity_
"""

# ╔═╡ e296c6e8-259c-11eb-1385-53f757f4d585
md"""
👉 Change the value of $B$ using the slider above. What does it mean for a climate system to have a more negative value of $B$? Explain why we call $B$ the _climate feedback parameter_.
"""

# ╔═╡ a86f13de-259d-11eb-3f46-1f6fb40020ce
observations_from_changing_B = md"""
Hello world!
"""

# ╔═╡ 3d66bd30-259d-11eb-2694-471fb3a4a7be
md"""
👉 What happens when $B$ is greater than or equal to zero?
"""

# ╔═╡ 5f82dec8-259e-11eb-2f4f-4d661f44ef41
observations_from_nonnegative_B = md"""
Hello world!
"""

# ╔═╡ 56b68356-2601-11eb-39a9-5f4b8e580b87
md"Reveal answer: $(@bind reveal_nonnegative_B_answer CheckBox())"

# ╔═╡ 7d815988-1fc7-11eb-322a-4509e7128ce3
if reveal_nonnegative_B_answer
	md"""
This is known as the "runaway greenhouse effect", where warming self-amplifies so strongly through *positive feedbacks* that the warming continues forever (or until the oceans boil away and there is no longer a reservoir or water to support a *water vapor feedback*. This is thought to explain Venus' extremely hot and hostile climate, but as you can see is extremely unlikely to occur on present-day Earth.
"""
end

# ╔═╡ 269200ec-259f-11eb-353b-0b73523ef71a
md"""
#### Exercise 1.2 - _Doubling CO₂_

To compute ECS, we doubled the CO₂ in our atmosphere. This factor 2 is not entirely arbitrary: without substantial effort to reduce CO₂ emissions, we are expected to **at least** double the CO₂ in our atmosphere by 2100. 

Right now, our CO₂ concentration is 415 ppm -- $(round(415 / 280, digits=3)) times the pre-industrial value of 280 ppm from 1850. 

The CO₂ concentrations in the _future_ depend on human action. There are several models for future emissions, which are formed by assuming different _policy scenarios_. A baseline model is RCP8.5 - a "worst-case" high-emissions scenario. In our notebook, this model is given as a function of ``t``.
"""

# ╔═╡ 2dfab366-25a1-11eb-15c9-b3dd9cd6b96c
md"""
👉 In what year are we expected to have doubled the CO₂ concentration, under policy scenario RCP8.5?
"""

# ╔═╡ 50ea30ba-25a1-11eb-05d8-b3d579f85652
expected_double_CO2_year = let
	
	
	missing
end

# ╔═╡ bade1372-25a1-11eb-35f4-4b43d4e8d156
md"""
The climate feedback parameter B is not something that we can control– it is an emergent property of the global climate system. Unfortunately, B is also difficult to quantify empirically (the relevant processes are difficult or impossible to observe directly), so there remains uncertainty as to its exact value.
"""

# ╔═╡ 02232964-2603-11eb-2c4c-c7b7e5fed7d1
B̅ = -1.3; σ = 0.4

# ╔═╡ c4398f9c-1fc4-11eb-0bbb-37f066c6027d
ECS(; B=B̅, a=Model.a) = -a*log(2.)./B;

# ╔═╡ 25f92dec-1fc4-11eb-055d-f34deea81d0e
let
	double_CO2(t) = if t >= 0
		2*Model.CO2_PI
	else
		Model.CO2_PI
	end
	
	# the definition of A depends on B, so we recalculate:
	A = Model.S*(1. - Model.α)/4 + B_slider*Model.T0
	# create the model
	ebm_ECS = Model.EBM(14., -100., 1., double_CO2, A=A, B=B_slider);
	Model.run!(ebm_ECS, 300)
	
	ecs = ECS(B=B_slider)
	
	p = plot(
		size=(500,250), legend=:bottomright, 
		title="Transient response to instant doubling of CO₂", 
		ylabel="temperature change [°C]", xlabel="years after doubling",
		ylim=(-.5, (isfinite(ecs) && ecs < 4) ? 4 : 10),
	)
	
	plot!(p, [ebm_ECS.t[1], ebm_ECS.t[end]], ecs .* [1,1], 
		ls=:dash, color=:darkred, label="ECS")
	
	plot!(p, ebm_ECS.t, ebm_ECS.T .- ebm_ECS.T[1], 
		label="ΔT(t) = T(t) - T₀")
end |> as_svg

# ╔═╡ 736ed1b6-1fc2-11eb-359e-a1be0a188670
B_samples = let
	B_distribution = Normal(B̅, σ)
	Nsamples = 5000
	
	rand(B_distribution, Nsamples)
end

# ╔═╡ 49cb5174-1fc3-11eb-3670-c3868c9b0255
histogram(B_samples, size=(600, 250), label=nothing, xlabel="B [W/m²/K]", ylabel="samples")

# ╔═╡ f3abc83c-1fc7-11eb-1aa8-01ce67c8bdde
md"""##### Exercise 1.4 - _Non-linear uncertainty propagation in climate_

**Question:** Use Monte Carlo simulation to generate a probability distribution for the ECS based on the probability distribution function for $B$ above.
"""

# ╔═╡ b6d7a362-1fc8-11eb-03bc-89464b55c6fc
md"**Answer:**"

# ╔═╡ 1f148d9a-1fc8-11eb-158e-9d784e390b24
begin
	ECS_samples = ECS.(B=B_samples);
	histogram(ECS_samples, xlims=(0, 8), size=(500, 240))
end

# ╔═╡ 6392bf28-210f-11eb-0793-835be433c454
scatter(B_samples, ECS_samples, ylims=[0, 20])

# ╔═╡ cf8dca6c-1fc8-11eb-1f89-099e6ba53c22
md"Compare the ECS distribution to the $\text{ECS}(\overline{B})$ that corresponds to the mean value of the climate feedback parameter $\overline{B}$.

👉 How does $\overline{\text{ECS}(B)}$ compare to $\text{ECS}(\overline{B})$? What is the probability that $\text{ECS}(B)$ lies above $\text{ECS}(\overline{B})$?
"

# ╔═╡ d44daea2-252f-11eb-364f-377ae504dc04
ecs_of_mean = ECS(B=mean(B_samples))

# ╔═╡ e27b2cd4-252f-11eb-20ef-0354db6220c2
mean_of_ecs = mean(ECS.(B=B_samples))

# ╔═╡ f94e635e-252f-11eb-1a52-310b628bd9b2
sum(ECS_samples) do e
	e > ecs_of_mean
end / length(ECS_samples)

# ╔═╡ 23e24d88-2530-11eb-26ef-c5e4e8b4f276
sum(ECS_samples) do e
	e > mean_of_ecs
end / length(ECS_samples)

# ╔═╡ 440271b6-25e8-11eb-26ce-1b80aa176aca
md"👉 Does accounting for uncertainty in feedbacks make our expectation of global warming better (less implied warming) or worse (more implied warming)?"

# ╔═╡ cf276892-25e7-11eb-38f0-03f75c90dd9e
observations_from_the_order_of_averaging = md"""
Hello world!
"""

# ╔═╡ 9c32db5c-1fc9-11eb-029a-d5d554de1067
md"""#### Exercise 1.5 - _Application to policy relevant questions_

We talked about two _emissions scenarios_: RCP2.6 (strong mitigation - controlled CO2 concentrations) and RCP8.5 (no mitigation - high CO2 concentrations). These are given by the following functions:
"""

# ╔═╡ ee1be5dc-252b-11eb-0865-291aa823b9e9
t = 1850:2100

# ╔═╡ e10a9b70-25a0-11eb-2aed-17ed8221c208
plot(t, Model.CO2_RCP85.(t), 
	ylim=(0,1200), ylabel="CO2 concentration [ppm]")

# ╔═╡ 40f1e7d8-252d-11eb-0549-49ca4e806e16
@bind t_scenario_test Slider(t; show_value=true, default=1850)

# ╔═╡ 19957754-252d-11eb-1e0a-930b5208f5ac
Model.CO2_RCP26(t_scenario_test), Model.CO2_RCP85(t_scenario_test)

# ╔═╡ 06c5139e-252d-11eb-2645-8b324b24c405
md"""
We are interested in how the **uncertainty in our input** $B$ (the climate feedback paramter) *propagates* through our model to determine the **uncertainty in our output** $T(t)$, for a given emissions scenario. The goal of this exercise is to answer the following by using *Monte Carlo Simulation* for *uncertainty propagation*:

> 👉 What is the probability that we see more than 2°C of warming by 2100 under the low-emissions scenario RCP2.6? What about under the high-emissions scenario RCP8.5?

"""

# ╔═╡ f2e55166-25ff-11eb-0297-796e97c62b07


# ╔═╡ 101cda5e-252e-11eb-2555-e3e8852f470f
md"""

**If Correct Answer:** shows a plot of the ''cone of uncertainty'' using `plot(t, T_low, fillrange=T_high)`
"""

# ╔═╡ 1ea81214-1fca-11eb-2442-7b0b448b49d6
md"""
## **Exercise 2** - _How did Snowball Earth melt?_

In lecture 21 (see below), we discovered that increases in the brightness of the Sun are not sufficient to explain how Snowball Earth eventually melted.
"""

# ╔═╡ a0ef04b0-25e9-11eb-1110-cde93601f712
html"""
<iframe width="100%" height="300" src="https://www.youtube-nocookie.com/embed/Y68tnH0FIzc" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 3e310cf8-25ec-11eb-07da-cb4a2c71ae34
md"""
We talked about a second theory -- a large increase in CO₂ (by volcanoes) could have caused a strong enough greenhouse effect to melt the Snowball. If we imagine that the CO₂ then decreased (e.g. by getting sequestered by the now liquid ocean), we might be able to explain how we transitioned from a hostile Snowball Earth to today's habitable "Waterball" Earth.

In this exercise, you will estimate how much CO₂ would be needed to melt the Snowball and visualize a possible trajectory for Earth's climate over the past 700 million years by making an interactive *bifurcation diagram*.

In the lecture notebook, we have a bifurcation diagram of $S$ (solar insolation) vs $T$ (temperature). We increased $S$, watched our point move right in the diagram until we found the tipping point. This time we will do the same, but we vary the CO₂ concentration, and keep $S$ fixed at its present day value.
"""

# ╔═╡ 0f52e312-2537-11eb-289e-17dc04710c2d
let
	ebm = Model.EBM(-40.0, 0., 5., t -> 280)
	
	Model.run!(ebm, 500)
	
	ebm.T
end

# ╔═╡ f984e274-2536-11eb-0092-27bb91984530
S = Model.S

# ╔═╡ 68b2a560-2536-11eb-0cc4-27793b4d6a70
function add_cold_hot_areas!(p)
	
	left, right = xlims(p)
	
	plot!([left, right], [-60, -60], fillrange=[-10., -10.], fillalpha=0.3, c=:lightblue, label=nothing)
	annotate!(left+12, -19, text("completely\nfrozen", 10, :darkblue, :left))
	
	plot!([left, right], [10, 10], fillrange=[80., 80.], fillalpha=0.09, c=:red, lw=0., label=nothing)
	annotate!(left+12, 15, text("no ice", 10, :darkred, :left))
end

# ╔═╡ c3e1deca-2530-11eb-0cb7-c3cc3118f1f6
begin
	CO2min = 10
	CO2max = 1_000_000
	
	CO2vec = CO2min:1.:CO2max
	CO2vec = vcat(CO2vec, reverse(CO2vec))
	Tvec = zeros(size(CO2vec))

	# local T_restart = -100.
	# for (i, CO2) = enumerate(CO2vec)
	# 	ebm = Model.EBM(T_restart, 0., 5., (t) -> CO2);
	# 	# ebm.S = S
	# 	Model.run!(ebm, 400.)
	# 	T_restart = ebm.T[end]
	# 	Tvec[i] = deepcopy(T_restart)
	# end
	
	md"**Data structures for storing warm & cool branch climates**"
end

# ╔═╡ 9f369200-2530-11eb-114c-6bb0bc2882af
# let
# 	ebm
# 	co2Slider = @bind CO2 Slider(CO2min:2.:CO2max, default=Model.CO2_PI);
# 	md""" $(CO2min) W/m² $(co2Slider) $(CO2max) W/m²"""
# end

# ╔═╡ 3c7d33da-253d-11eb-0c5a-9b0d524c42f8
@bind log_CO2 Slider(log10(CO2min):0.01:log10(CO2max); default=log10(Model.CO2_PI))

# ╔═╡ 35f87c2e-253d-11eb-0d79-61d89c1d9b5e
CO2 = 10^log_CO2

# ╔═╡ aa1a3562-2537-11eb-0010-abde7b40090a
function restart_ebm!(ebm)
	ebm.T = [ebm.T[end]]
	ebm.t = [ebm.t[1]]
end

# ╔═╡ e411a3bc-2538-11eb-3492-bfdd42b1445d
function step_model!(ebm, new_CO2)
	restart_ebm!(ebm)
	
	
	ebm.CO2 = t -> new_CO2
	Model.run!(ebm, 500)
	
	ebm
end

# ╔═╡ d7801e88-2530-11eb-0b93-6f1c78d00eea
function α(T; α0=Model.α, αi=0.5, ΔT=10.)
	if T < -ΔT
		return αi
	elseif -ΔT <= T < ΔT
		return αi + (α0-αi)*(T+ΔT)/(2ΔT)
	elseif T >= ΔT
		return α0
	end
end

# ╔═╡ 607058ec-253c-11eb-0fb6-add8cfb73a4f
function Model.timestep!(ebm)
	ebm.α = α(ebm.T[end]) # Added this line
	append!(ebm.T, ebm.T[end] + ebm.Δt*Model.tendency(ebm));
	append!(ebm.t, ebm.t[end] + ebm.Δt);
end

# ╔═╡ cdc54b98-2530-11eb-3d5e-71c4b53256fb
begin
	Sneo = Model.S*0.93
	Tneo = -48.
	md"**Initial conditions**"
end

# ╔═╡ 06d28052-2531-11eb-39e2-e9613ab0401c
begin
	ebm = Model.EBM(Tneo, 0., 5., Model.CO2_const)
	
	md"**Data structures for storing trajectories of recent climates**"
end

# ╔═╡ fc94977a-253b-11eb-1143-912bb5ef055f
ebm.α

# ╔═╡ 378aed18-252b-11eb-0b37-a3b511af2cb5
let
	CO2
	
	step_model!(ebm, CO2)
	
	p = plot(
		xlims=(CO2min, CO2max), ylims=(-55, 75), 
		xaxis=:log,
		title="Earth's CO2 concentration bifurcation diagram"
	)
	plot!([Model.CO2_PI, Model.CO2_PI], [-55, 75], color=:grey, alpha=0.3, lw=8, label="Pre-industrial CO2")
	if false
		plot!(p, xlims=(CO2min, CO2max))
		# if show_cold
		# 	plot!(CO2vec[warming_mask], Tvec[warming_mask], color=:blue,lw=3., alpha=0.5, label="cool branch")
		# end
		# if show_warm
		# 	plot!(CO2vec[.!warming_mask], Tvec[.!warming_mask], color="red", lw=3., alpha=0.5, label="warm branch")
		# end
		# if show_unstable
		# 	plot!(CO_unstable, T_unstable, color=:darkgray, lw=3., alpha=0.4, ls=:dash, label="unstable branch")
		# end
	end
	plot!(legend=:topleft)
	plot!(xlabel="CO2 concentration [ppm]", ylabel="Global temperature T [°C]")
	plot!([Model.CO2_PI], [Model.T0], shape=:circle, label="Our preindustrial climate", color=:orange, markersize=8)
	plot!([Model.CO2_PI], [-38.3], shape=:circle, label="Alternate preindustrial climate", color=:aqua, markersize=8)
	# plot!([Sneo], [Tneo], marker=:., label="neoproterozoic (700 Mya)", color=:lightblue, markersize=8)
	z = [ebm.CO2(123), ebm.T[end]]
	# plot!(plot!(), z, color="black", marker=:c, markersize=8/2*1.2, label=nothing, markerstrokecolor=nothing, markerstrokewidth=0.)
	scatter!(z[1:1], z[2:2])
	
	add_cold_hot_areas!(plot!())
end |> as_svg

# ╔═╡ 2f39805a-25ed-11eb-0594-89755c16ad15
CO2vec_hires = CO2min:0.1:CO2max

# ╔═╡ cb15cd88-25ed-11eb-2be4-f31500a726c8
md"Hint: Use a condition on the albedo or temperature to check whether the Snowball has melted."

# ╔═╡ 232b9bec-2544-11eb-0401-97a60bb172fc
md"Hint: Start by writing a function `equilibrate(CO2)` which starts at the Snowball Earth temperature T = $(Tneo) and returns the equilibrium temperature for a given CO2 level."

# ╔═╡ 1dcce868-2544-11eb-14af-4f7811b7f2a8


# ╔═╡ 3a35598a-2527-11eb-37e5-3b3e4c63c4f7
md"""
## **Exercise XX:** _Lecture transcript_
_(MIT students only)_

Please see the link for hw 9 transcript document on [Canvas](https://canvas.mit.edu/courses/5637).
We want each of you to correct about 500 lines, but don’t spend more than 20 minutes on it.
See the the beginning of the document for more instructions.
:point_right: Please mention the name of the video(s) and the line ranges you edited:
"""

# ╔═╡ 5041cdee-2527-11eb-154f-0b0c68e11fe3
lines_i_edited = md"""
Abstraction, lines 1-219; Array Basics, lines 1-137; Course Intro, lines 1-144 (_for example_)
"""

# ╔═╡ 36e2dfea-2433-11eb-1c90-bb93ab25b33c
if student.name == "Jazzy Doe" || student.kerberos_id == "jazz"
	md"""
	!!! danger "Before you submit"
	    Remember to fill in your **name** and **Kerberos ID** at the top of this notebook.
	"""
end

# ╔═╡ 36ea4410-2433-11eb-1d98-ab4016245d95
md"## Function library

Just some helper functions used in the notebook."

# ╔═╡ 36f8c1e8-2433-11eb-1f6e-69dc552a4a07
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))

# ╔═╡ 51e2e742-25a1-11eb-2511-ab3434eacc3e
hint(md"The function `findfirst` might be helpful.")

# ╔═╡ 37061f1e-2433-11eb-3879-2d31dc70a771
almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]))

# ╔═╡ 371352ec-2433-11eb-153d-379afa8ed15e
still_missing(text=md"Replace `missing` with your answer.") = Markdown.MD(Markdown.Admonition("warning", "Here we go!", [text]))

# ╔═╡ 372002e4-2433-11eb-0b25-39ce1b1dd3d1
keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]))

# ╔═╡ 372c1480-2433-11eb-3c4e-95a37d51835f
yays = [md"Fantastic!", md"Splendid!", md"Great!", md"Yay ❤", md"Great! 🎉", md"Well done!", md"Keep it up!", md"Good job!", md"Awesome!", md"You got the right answer!", md"Let's move on to the next section."]

# ╔═╡ 3737be8e-2433-11eb-2049-2d6d8a5e4753
correct(text=rand(yays)) = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]))

# ╔═╡ 374522c4-2433-11eb-3da3-17419949defc
not_defined(variable_name) = Markdown.MD(Markdown.Admonition("danger", "Oopsie!", [md"Make sure that you define a variable called **$(Markdown.Code(string(variable_name)))**"]))

# ╔═╡ 37552044-2433-11eb-1984-d16e355a7c10
TODO = html"<span style='display: inline; font-size: 2em; color: purple; font-weight: 900;'>TODO</span>"

# ╔═╡ 291326e8-25a2-11eb-1a00-3de0f60e5f0f
md"""
#### Exercise 1.3 - _Uncertainty in B_

$TODO

The point of this exercise is:

1. there is a small chance that B is close to zero
1. but B close to zero has dramatic effects
1. 👉 nonlinearity on distributions
"""

# ╔═╡ a2aff256-1fc6-11eb-3671-b7801bce27fc
md"""In Exercise 1 we talked about the meaning of ``B \geq 0``. How likely is this scenario?

$TODO does the original paper say anything about this B>0 tail of the distribution? They might not have intended to assign a probability to B>0

No, they always clip it at zero but they'd explain why. I think it's still worth talking about it because it brings up the concept of a "runaway feedback". I think we should still include this in the problem set but we can clarify that the tails of the climate sensivity distribution are less well constrained than the center– i.e. it's possible that the actual distribution is non-Gaussian.
"""

# ╔═╡ d6d1b312-2543-11eb-1cb2-e5b801686ffb
md"""
$TODO

the exercise here will be that you start out with an "empty" CO2 plot, and you add the ebm with CO2 slider.
"""

# ╔═╡ f81da5ee-2543-11eb-0f34-93b47dbf4c34
md"""
$TODO

exercise: add the trail to this viz using push! and pop! (or using a circular buffer)
"""

# ╔═╡ 11096250-2544-11eb-057b-d7112f20b05c
md"""
$TODO

Find the **lowest CO₂ concentration** necessary to melt the Snowball, programatically, by performing a *binary search* on `CO2vec_hires`?
"""

# ╔═╡ Cell order:
# ╟─169727be-2433-11eb-07ae-ab7976b5be90
# ╟─18be4f7c-2433-11eb-33cb-8d90ca6f124c
# ╟─21524c08-2433-11eb-0c55-47b1bdc9e459
# ╠═23335418-2433-11eb-05e4-2b35dc6cca0e
# ╟─253f4da0-2433-11eb-1e48-4906059607d3
# ╠═1e06178a-1fbf-11eb-32b3-61769a79b7c0
# ╟─87e68a4a-2433-11eb-3e9d-21675850ed71
# ╟─1312525c-1fc0-11eb-2756-5bc3101d2260
# ╠═c4398f9c-1fc4-11eb-0bbb-37f066c6027d
# ╟─7f961bc0-1fc5-11eb-1f18-612aeff0d8df
# ╟─930d7154-1fbf-11eb-1c3a-b1970d291811
# ╟─25f92dec-1fc4-11eb-055d-f34deea81d0e
# ╟─fa7e6f7e-2434-11eb-1e61-1b1858bb0988
# ╟─16348b6a-1fc2-11eb-0b9c-65df528db2a1
# ╟─e296c6e8-259c-11eb-1385-53f757f4d585
# ╠═a86f13de-259d-11eb-3f46-1f6fb40020ce
# ╟─3d66bd30-259d-11eb-2694-471fb3a4a7be
# ╠═5f82dec8-259e-11eb-2f4f-4d661f44ef41
# ╟─56b68356-2601-11eb-39a9-5f4b8e580b87
# ╟─7d815988-1fc7-11eb-322a-4509e7128ce3
# ╟─269200ec-259f-11eb-353b-0b73523ef71a
# ╠═e10a9b70-25a0-11eb-2aed-17ed8221c208
# ╟─2dfab366-25a1-11eb-15c9-b3dd9cd6b96c
# ╠═50ea30ba-25a1-11eb-05d8-b3d579f85652
# ╟─51e2e742-25a1-11eb-2511-ab3434eacc3e
# ╟─291326e8-25a2-11eb-1a00-3de0f60e5f0f
# ╟─bade1372-25a1-11eb-35f4-4b43d4e8d156
# ╠═02232964-2603-11eb-2c4c-c7b7e5fed7d1
# ╠═736ed1b6-1fc2-11eb-359e-a1be0a188670
# ╠═49cb5174-1fc3-11eb-3670-c3868c9b0255
# ╠═a2aff256-1fc6-11eb-3671-b7801bce27fc
# ╠═6392bf28-210f-11eb-0793-835be433c454
# ╟─f3abc83c-1fc7-11eb-1aa8-01ce67c8bdde
# ╟─b6d7a362-1fc8-11eb-03bc-89464b55c6fc
# ╠═1f148d9a-1fc8-11eb-158e-9d784e390b24
# ╠═cf8dca6c-1fc8-11eb-1f89-099e6ba53c22
# ╠═d44daea2-252f-11eb-364f-377ae504dc04
# ╠═e27b2cd4-252f-11eb-20ef-0354db6220c2
# ╠═f94e635e-252f-11eb-1a52-310b628bd9b2
# ╠═23e24d88-2530-11eb-26ef-c5e4e8b4f276
# ╟─440271b6-25e8-11eb-26ce-1b80aa176aca
# ╠═cf276892-25e7-11eb-38f0-03f75c90dd9e
# ╟─9c32db5c-1fc9-11eb-029a-d5d554de1067
# ╠═19957754-252d-11eb-1e0a-930b5208f5ac
# ╠═40f1e7d8-252d-11eb-0549-49ca4e806e16
# ╟─ee1be5dc-252b-11eb-0865-291aa823b9e9
# ╟─06c5139e-252d-11eb-2645-8b324b24c405
# ╠═f2e55166-25ff-11eb-0297-796e97c62b07
# ╠═101cda5e-252e-11eb-2555-e3e8852f470f
# ╟─1ea81214-1fca-11eb-2442-7b0b448b49d6
# ╟─a0ef04b0-25e9-11eb-1110-cde93601f712
# ╟─3e310cf8-25ec-11eb-07da-cb4a2c71ae34
# ╠═d6d1b312-2543-11eb-1cb2-e5b801686ffb
# ╠═0f52e312-2537-11eb-289e-17dc04710c2d
# ╠═fc94977a-253b-11eb-1143-912bb5ef055f
# ╠═f984e274-2536-11eb-0092-27bb91984530
# ╠═68b2a560-2536-11eb-0cc4-27793b4d6a70
# ╠═c3e1deca-2530-11eb-0cb7-c3cc3118f1f6
# ╠═9f369200-2530-11eb-114c-6bb0bc2882af
# ╠═3c7d33da-253d-11eb-0c5a-9b0d524c42f8
# ╠═35f87c2e-253d-11eb-0d79-61d89c1d9b5e
# ╠═e411a3bc-2538-11eb-3492-bfdd42b1445d
# ╠═378aed18-252b-11eb-0b37-a3b511af2cb5
# ╠═06d28052-2531-11eb-39e2-e9613ab0401c
# ╠═aa1a3562-2537-11eb-0010-abde7b40090a
# ╠═d7801e88-2530-11eb-0b93-6f1c78d00eea
# ╠═607058ec-253c-11eb-0fb6-add8cfb73a4f
# ╠═cdc54b98-2530-11eb-3d5e-71c4b53256fb
# ╠═f81da5ee-2543-11eb-0f34-93b47dbf4c34
# ╠═11096250-2544-11eb-057b-d7112f20b05c
# ╠═2f39805a-25ed-11eb-0594-89755c16ad15
# ╟─cb15cd88-25ed-11eb-2be4-f31500a726c8
# ╠═232b9bec-2544-11eb-0401-97a60bb172fc
# ╠═1dcce868-2544-11eb-14af-4f7811b7f2a8
# ╟─3a35598a-2527-11eb-37e5-3b3e4c63c4f7
# ╠═5041cdee-2527-11eb-154f-0b0c68e11fe3
# ╟─36e2dfea-2433-11eb-1c90-bb93ab25b33c
# ╟─36ea4410-2433-11eb-1d98-ab4016245d95
# ╟─36f8c1e8-2433-11eb-1f6e-69dc552a4a07
# ╟─37061f1e-2433-11eb-3879-2d31dc70a771
# ╟─371352ec-2433-11eb-153d-379afa8ed15e
# ╟─372002e4-2433-11eb-0b25-39ce1b1dd3d1
# ╟─372c1480-2433-11eb-3c4e-95a37d51835f
# ╟─3737be8e-2433-11eb-2049-2d6d8a5e4753
# ╟─374522c4-2433-11eb-3da3-17419949defc
# ╟─37552044-2433-11eb-1984-d16e355a7c10
