### A Pluto.jl notebook ###
# v0.12.9

using Markdown
using InteractiveUtils

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
md"""## Problem 1: policy goals under uncertainty
A recent ground-breaking review paper produced the most comprehensive and up-to-date estimate of the *climate feedback parameter*, which they find to be

$B \approx \mathcal{N}(-1.3, 0.4),$

i.e. is normally distributed with a mean value $\overline{B} = -1.3$ W/m²/K and a standard deviation $\sigma = 0.4$ W/m²/K. These value are not very intuitive, so let us convert them into more policy-relevant numbers.

**Definition:** *Equilibrium climate sensitivity (ECS)* is defined as the amount of warming $\Delta T$ caused by a doubling of CO₂ (e.g. from the pre-industrial value 280 ppm to 560 ppm), at equilibrium.

At equilibrium, the energy balance model equation is:

$0 = \frac{S(1 - α)}{4} - (A - BT_{eq}) + a \ln\left( \frac{2\;\text{CO}₂_{\text{PI}}}{\text{CO}₂_{\text{PI}}} \right)$

Subtracting the preindustrial energy balance 

$0 = \frac{S(1-α)}{4} - (A - BT_{0}),$

we have

$\text{ECS} \equiv T_{eq} - T_{0} = -\frac{a\ln(2)}{B}$
"""

# ╔═╡ 7f961bc0-1fc5-11eb-1f18-612aeff0d8df
md"""The plot below provides an example of an "abrupt 2xCO₂" experiment, a classic experimental treatment method in climate modelling which is used in practice to estimate ECS for a particular model (Note: in complicated climate models the values of the parameters $a$ and $B$ are not specified *apriori*, but emerge as outputs for the simulation).

The simulation begins at the preindustrial equilibrium, i.e. a temperature $T_{0} = 14$°C is in balance with the pre-industrial CO₂ concentration of 280 ppm until CO₂ is abruptly doubled from 280 ppm to 560 ppm. The climate responds by rapidly warming, and after a few hundred years approaches the equilibrium climate sensitivity value, by definition.
"""

# ╔═╡ 16348b6a-1fc2-11eb-0b9c-65df528db2a1
md"""
##### Problem 1. (a) Develop understanding for feedbacks and climate sensitivity
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

# ╔═╡ 736ed1b6-1fc2-11eb-359e-a1be0a188670
begin
	B̅ = -1.3; σ = 0.4
	d = Normal(B̅, σ)
	Nsamples = 5000
	
	B_samples = rand(d, Nsamples)
end;

# ╔═╡ c4398f9c-1fc4-11eb-0bbb-37f066c6027d
ECS(; B=B̅, a=Model.a) = -a*log(2.)./B;

# ╔═╡ 25f92dec-1fc4-11eb-055d-f34deea81d0e
begin
	double_CO2(t) = 2*Model.CO2_PI
	ebm_ECS = Model.EBM(14., 0., 1., double_CO2, B=B̅);
	Model.run!(ebm_ECS, 300)
	plot(size=(500,250), legend=:bottomright, title="Transient response to instant doubling of CO₂", ylabel="temperature [°C]", xlabel="years after doubling")
	plot!([0, 300], ECS() .* [1,1], ls=:dash, color=:darkred, label="ECS")
	plot!(ebm_ECS.t, ebm_ECS.T .- ebm_ECS.T[1], label="ΔT(t) = T(t) - T₀")
end

# ╔═╡ 49cb5174-1fc3-11eb-3670-c3868c9b0255
histogram(B_samples, size=(600, 250), label=nothing, xlabel="B [W/m²/K]", ylabel="samples")

# ╔═╡ a2aff256-1fc6-11eb-3671-b7801bce27fc
md"**Question:** What happens if the climate feedback parameter $B$ is less than or equal to zero? How often does this scenario occur?"

# ╔═╡ 82f8fe38-1fc3-11eb-3a89-ffe737246a28
begin
	ebm = Model.EBM(14., 0., 1., Model.CO2_const, B=0.);
	Model.run!(ebm, 500)
	plot(ebm.t, ebm.T, size=(300, 250), ylabel="temperature [°C]", xlabel="year", label=nothing)
end

# ╔═╡ 7d815988-1fc7-11eb-322a-4509e7128ce3
md"""**Answer:** endless warming!!! ahhhhh

**If answered correctly:** This is known as the "runaway greenhouse effect", where warming self-amplifies so strongly through *positive feedbacks* that the warming continues forever (or until the oceans boil away and there is no longer a reservoir or water to support a *water vapor feedback*. This is thought to explain Venus' extremely hot and hostile climate, but as you can see is extremely unlikely to occur on present-day Earth.
"""

# ╔═╡ f3abc83c-1fc7-11eb-1aa8-01ce67c8bdde
md"""##### Problem 2. (b) Non-linear uncertainty propagation in climate

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
md"**Question:** Compare the ECS distribution to the $\text{ECS}(\overline{B})$ that corresponds to the mean value of the climate feedback parameter $\overline{B}$.

How does $\overline{\text{ECS}(B)}$ compare to $\text{ECS}(\overline{B})$? What is the probability that $\text{ECS}(B)$ lies above $\text{ECS}(\overline{B})$?"

# ╔═╡ 9c32db5c-1fc9-11eb-029a-d5d554de1067
md"##### Problem 1. (c) Application to policy relevant questions

**Question:** What is the probability that we see more than 2°C of warming by 2100 under the low-emissions scenario RCP2.6? What about under the high-emissions scenario RCP8.5?

**If Correct Answer:** shows a plot of the ''cone of uncertainty'' using `plot(t, T_low, fillrange=T_high)`
"

# ╔═╡ 1ea81214-1fca-11eb-2442-7b0b448b49d6
md"""
## Problem 2. How did Snowball Earth melt?

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
# ╠═25f92dec-1fc4-11eb-055d-f34deea81d0e
# ╟─16348b6a-1fc2-11eb-0b9c-65df528db2a1
# ╟─930d7154-1fbf-11eb-1c3a-b1970d291811
# ╠═736ed1b6-1fc2-11eb-359e-a1be0a188670
# ╟─49cb5174-1fc3-11eb-3670-c3868c9b0255
# ╟─a2aff256-1fc6-11eb-3671-b7801bce27fc
# ╠═82f8fe38-1fc3-11eb-3a89-ffe737246a28
# ╠═6392bf28-210f-11eb-0793-835be433c454
# ╟─7d815988-1fc7-11eb-322a-4509e7128ce3
# ╟─f3abc83c-1fc7-11eb-1aa8-01ce67c8bdde
# ╟─b6d7a362-1fc8-11eb-03bc-89464b55c6fc
# ╠═1f148d9a-1fc8-11eb-158e-9d784e390b24
# ╟─cf8dca6c-1fc8-11eb-1f89-099e6ba53c22
# ╟─9c32db5c-1fc9-11eb-029a-d5d554de1067
# ╟─1ea81214-1fca-11eb-2442-7b0b448b49d6
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
