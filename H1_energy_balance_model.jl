### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ 1e06178a-1fbf-11eb-32b3-61769a79b7c0
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	Pkg.add("LaTeXStrings")
	Pkg.add("Distributions")
	Pkg.add("Random")
	using LaTeXStrings
	using Plots
	using PlutoUI
	using Random, Distributions
	Random.seed!(123)
	
	md"##### Package dependencies"
end

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

# ╔═╡ 736ed1b6-1fc2-11eb-359e-a1be0a188670
begin
	B̅ = -1.3; σ = 0.4
	d = Normal(B̅, σ)
	Nsamples = 5000
	
	B_samples = rand(d, Nsamples)
end;

# ╔═╡ 49cb5174-1fc3-11eb-3670-c3868c9b0255
histogram(B_samples, size=(600, 250), label=nothing, xlabel="B [W/m²/K]", ylabel="samples")

# ╔═╡ a2aff256-1fc6-11eb-3671-b7801bce27fc
md"**Question:** What happens if the climate feedback parameter $B$ is less than or equal to zero? How often does this scenario occur?"

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

# ╔═╡ 616a0bda-1fbf-11eb-3263-216d5853f8a5
md"""## Pluto magic"""

# ╔═╡ 3ebd3ab2-1fbf-11eb-0a7f-05c82b8013c4
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 930d7154-1fbf-11eb-1c3a-b1970d291811
Model = ingredients("1_energy_balance_model.jl");

# ╔═╡ c4398f9c-1fc4-11eb-0bbb-37f066c6027d
ECS(; B=B̅, a=Model.a) = -a*log(2.)./B;

# ╔═╡ 1f148d9a-1fc8-11eb-158e-9d784e390b24
begin
	ECS_samples = ECS.(B=B_samples);
	histogram(ECS_samples, xlims=(0, 8), size=(500, 240))
end

# ╔═╡ 6392bf28-210f-11eb-0793-835be433c454
scatter(B_samples, ECS_samples, ylims=[0, 20])

# ╔═╡ 25f92dec-1fc4-11eb-055d-f34deea81d0e
begin
	double_CO2(t) = 2*Model.CO2_PI
	ebm_ECS = Model.EBM(14., 0., 1., double_CO2, B=B̅);
	Model.run!(ebm_ECS, 300)
	plot(size=(500,250), legend=:bottomright, title="Transient response to instant doubling of CO₂", ylabel="temperature [°C]", xlabel="years after doubling")
	plot!([0, 300], ECS() .* [1,1], ls=:dash, color=:darkred, label="ECS")
	plot!(ebm_ECS.t, ebm_ECS.T .- ebm_ECS.T[1], label="ΔT(t) = T(t) - T₀")
end

# ╔═╡ 82f8fe38-1fc3-11eb-3a89-ffe737246a28
begin
	ebm = Model.EBM(14., 0., 1., Model.CO2_const, B=0.);
	Model.run!(ebm, 500)
	plot(ebm.t, ebm.T, size=(300, 250), ylabel="temperature [°C]", xlabel="year", label=nothing)
end

# ╔═╡ Cell order:
# ╟─1312525c-1fc0-11eb-2756-5bc3101d2260
# ╠═c4398f9c-1fc4-11eb-0bbb-37f066c6027d
# ╟─7f961bc0-1fc5-11eb-1f18-612aeff0d8df
# ╠═25f92dec-1fc4-11eb-055d-f34deea81d0e
# ╟─16348b6a-1fc2-11eb-0b9c-65df528db2a1
# ╠═930d7154-1fbf-11eb-1c3a-b1970d291811
# ╠═736ed1b6-1fc2-11eb-359e-a1be0a188670
# ╟─49cb5174-1fc3-11eb-3670-c3868c9b0255
# ╟─a2aff256-1fc6-11eb-3671-b7801bce27fc
# ╟─82f8fe38-1fc3-11eb-3a89-ffe737246a28
# ╠═6392bf28-210f-11eb-0793-835be433c454
# ╟─7d815988-1fc7-11eb-322a-4509e7128ce3
# ╟─f3abc83c-1fc7-11eb-1aa8-01ce67c8bdde
# ╟─b6d7a362-1fc8-11eb-03bc-89464b55c6fc
# ╠═1f148d9a-1fc8-11eb-158e-9d784e390b24
# ╟─cf8dca6c-1fc8-11eb-1f89-099e6ba53c22
# ╟─9c32db5c-1fc9-11eb-029a-d5d554de1067
# ╟─1ea81214-1fca-11eb-2442-7b0b448b49d6
# ╟─616a0bda-1fbf-11eb-3263-216d5853f8a5
# ╟─3ebd3ab2-1fbf-11eb-0a7f-05c82b8013c4
# ╠═1e06178a-1fbf-11eb-32b3-61769a79b7c0
