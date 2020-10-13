### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ e9ad66b0-0d6b-11eb-26c0-27413c19dd32
begin
	using Pkg
	Pkg.add("Plots")
	using Plots
end

# ╔═╡ fbf1ba3e-0d65-11eb-20a7-55402d46d4ed
md"""
### A simple energy balance model

$C \frac{dT}{dt} = \frac{(1 - α)S}{4} - (A + BT).$

Things become much more interesting when we consider the ice-albedo feedback, i.e. the fact that the reflectively (or *albedo*) of the planet depends strongly on whether it is mostly covered in liquid water (which is a very absorbent dark blue) or solid ice (which is a very reflective white). We can approximate this behavior globally by specifying the following step-function for the albedo's dependence on temperature:

$\alpha(T) = \begin{cases}
\alpha_{o} &\mbox{if } T > 0\text{°C  (water in liquid phase)}\\
\alpha_{i} & \mbox{if } T \leq 0\text{°C  (water in solid phase)}
\end{cases}$

which can be **discretized** in time as

$C \frac{T_{n+1} - T_{n}}{\Delta t} = \frac{1-\alpha(T_{n}) S}{4} - (A + BT_{n}),$

where $n$ is the present timestep and $n+1$ is the next timestep $\Delta t$ later, which we can rewrite as

$T_{n+1} = \Delta t \left[ \frac{1}{C} \left( \frac{1-\alpha(T_{n}) S}{4} - (A + BT_{n}) \right) \right]$
"""

# ╔═╡ 6b57e7f8-0d7b-11eb-272e-010167db303b
md"""
We begin by creating a data structure to store the present temperature and timestep.
"""

# ╔═╡ a3d177b4-0d67-11eb-3578-d1f2d94f10b2
begin
	mutable struct EBM
		T::Float64
		t::Float64
	end
	EBM(T) = EBM(T, 0.)
end

# ╔═╡ 7e151456-0d7b-11eb-3981-39f1dc096a49
md"""Then, we set the values of the various parameters"""

# ╔═╡ 47fce926-0d66-11eb-1659-a1d6e88b3251

begin
	dt = 1.;
	
	C = 40.
	αo = 0.30
	αi = 0.20
	S0 = 1365.
	
	T0 = 14.
	B = 1.13
	A = (1. - αo) * S0/4 - B*T0
end;

# ╔═╡ a20699f2-0d7b-11eb-3e35-3923eb70affe
md"Next, we specify the function that determines the albedo for a given temperature"

# ╔═╡ 9dcae4ba-0d7b-11eb-296e-9745d3f5d70d
α(T) = αi*(T < 0.) + αo;

# ╔═╡ bdc26c98-0d7b-11eb-039c-59253d5f766d
md"And finally we define our `timestep!` method, which is what we use to solve the energy balance model."

# ╔═╡ 283abb98-0d66-11eb-14a0-f9673b11f38d
function timestep!(ebm::EBM; S=S0)
	ebm.T += dt * (1. /C) * ( (1 - α(ebm.T))*S/4 - (A + B*ebm.T))
	ebm.t += dt
end

# ╔═╡ fdabcf98-0d7b-11eb-0ca8-39787e1eedd9
md"As a shortcut, we'll also define a `run!` method that will run our energy balance model for a specified amount of time."

# ╔═╡ dc68c934-0d68-11eb-02cc-5d181cfaaae3
function run!(ebm; nt=1000, S=S0)
	for i=1.:nt
		timestep!(ebm, S=S)
	end
end

# ╔═╡ 29f90d22-0d7c-11eb-3a04-87ccea95d4f1
md"""
#### Exploring Earth's Multiple Equilibria

About 700 million years ago, the Earth underwent several dramatic climate changes and rapidly plunged from a warm, wet climate into a "snowball"-like state wherein the ocean froze all the way from the poles to the tropics.

How could this have happened?

"""

# ╔═╡ b5849000-0d68-11eb-3beb-c575e8d0ce8e
begin
	ebm = EBM(-100.)
	
	Svec = 1200.:1.:1850.
	Svec = vcat(Svec, reverse(Svec))
	Tvec = zeros(size(Svec))
	
	for (i, S)=enumerate(Svec)
		run!(ebm, S=S)
		Tvec[i] = copy(ebm.T)
	end
end;

# ╔═╡ 0f222836-0d6c-11eb-2ee8-45545da73cfd
begin
	warming_mask = (1:size(Svec)[1]) .< size(Svec)./2
	plot(Svec[warming_mask], Tvec[warming_mask], color="blue",lw=3., alpha=0.5, label="cool branch")
	plot!(Svec[.!warming_mask], Tvec[.!warming_mask], color="red", lw=3., alpha=0.5, label="warm branch")
	plot!(legend=:topleft, xlabel="solar insolation S₀ [W/m²]", ylabel="Global temperature [°C]")
	plot!([S0], [T0], marker=:c, label="present day", color="black", markersize=8)
	plot!([S0*0.92], [-58.3], marker=:c, label="neoproterozoic", color="lightblue", markersize=8)
end

# ╔═╡ 5cf4ccdc-0d7e-11eb-2683-fd9a72e763f2
md"""
### Exercises:
1. What happens if we instead assume the albedo varies continusouly from -10°C to 10°C? Repeat the above analysis with

$\alpha(T) = \begin{cases}
\alpha_{o} & \mbox{if } T < -10°\text{C}\\
\alpha_{o} + \alpha_{i} \frac{(T+10)}{20} & \mbox{if } -10°\text{C}<T<10°\text{C} \\
\alpha_{i} & \mbox{if } T > 10°\text{C}\\
\end{cases}$

2. What happens if we add CO2 to the neoprotorezoic atmosphere? How much CO2 would we need to add (relevant to a modern background of 280 ppm) to un-freeze the snowball?
"""

# ╔═╡ Cell order:
# ╟─fbf1ba3e-0d65-11eb-20a7-55402d46d4ed
# ╟─6b57e7f8-0d7b-11eb-272e-010167db303b
# ╠═a3d177b4-0d67-11eb-3578-d1f2d94f10b2
# ╟─7e151456-0d7b-11eb-3981-39f1dc096a49
# ╠═47fce926-0d66-11eb-1659-a1d6e88b3251
# ╟─a20699f2-0d7b-11eb-3e35-3923eb70affe
# ╠═9dcae4ba-0d7b-11eb-296e-9745d3f5d70d
# ╠═bdc26c98-0d7b-11eb-039c-59253d5f766d
# ╠═283abb98-0d66-11eb-14a0-f9673b11f38d
# ╟─fdabcf98-0d7b-11eb-0ca8-39787e1eedd9
# ╠═dc68c934-0d68-11eb-02cc-5d181cfaaae3
# ╟─29f90d22-0d7c-11eb-3a04-87ccea95d4f1
# ╠═b5849000-0d68-11eb-3beb-c575e8d0ce8e
# ╠═0f222836-0d6c-11eb-2ee8-45545da73cfd
# ╟─5cf4ccdc-0d7e-11eb-2683-fd9a72e763f2
# ╠═e9ad66b0-0d6b-11eb-26c0-27413c19dd32
