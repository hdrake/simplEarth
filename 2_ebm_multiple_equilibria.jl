### A Pluto.jl notebook ###
# v0.12.6

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

# ╔═╡ 05031b60-1df4-11eb-2b61-956e526b3d4a
md"## Lecture 2: Snowball Earth, the ice-albedo feedback, and multiple equilibria"

# ╔═╡ 016c1074-1df4-11eb-2da8-578e25d9456b
md"### The non-linear ice-albedo feedback"

# ╔═╡ 9c118f9a-1df0-11eb-22dd-b14428994076
md"![](https://upload.wikimedia.org/wikipedia/commons/d/df/Ice_albedo_feedback.jpg)"

# ╔═╡ 38346e6a-0d98-11eb-280b-f79787a3c788
md"""

For this section, we will ignore the influence of humans and thus assume the CO₂ concentrations remain equal to their pre-industrial values.

The *ice-albedo feedback* refers to the fact that the reflectively (or *albedo*) of the planet depends strongly on whether it is mostly covered in liquid water (which is a very absorbent dark blue) or solid ice (which is a very reflective white). We can approximate this behavior globally by specifying the following step-function for the albedo's dependence on temperature:

$\alpha(T) = \begin{cases}
\alpha_{0} &\mbox{if } T > 0\text{°C  (water in liquid phase)}\\
\alpha_{i} & \mbox{if } T \leq 0\text{°C  (water in solid phase)}
\end{cases}$
"""

# ╔═╡ a8dcc0fc-1df8-11eb-21fd-1fdebe5dabfc
md"""
All we have to do to modify our energy balance model from Lecture 1 is to overwrite the definition of the timestep! to specify that the albedo be a function of the present temperature.
"""

# ╔═╡ 29f90d22-0d7c-11eb-3a04-87ccea95d4f1
md"""
#### Exploring Earth's Multiple Equilibria

About 700 million years ago, the Earth underwent several dramatic climate changes and rapidly plunged from a warm, wet climate into a "snowball"-like state wherein the ocean froze all the way from the poles to the tropics.

How could this have happened?

"""

# ╔═╡ 61fc91ec-1df8-11eb-13c1-098c113b46ec
function restart_ebm!(ebm)
	ebm.T = [ebm.T[end]]
	ebm.t = [ebm.t[1]]
end

# ╔═╡ 9b39df12-1df9-11eb-2eb0-f138980be597
function plot_trajectory!(p, x, y; lw=8)
	n = size(x,1)
	plot!(x, y, color="black", alpha=collect(1/n:1/n:1.), linewidth=collect(0.:lw/n:lw), label=nothing)
	plot!((x[end], y[end]), color="black", marker=:c, markersize=lw/2*1.2, label=nothing, markerstrokecolor=nothing, markerstrokewidth=0.)
	return p
end;

# ╔═╡ 90ae18dc-0db8-11eb-0d73-c3b7efaef9b0
begin
		Sneo = 1250.
		Tneo = -50.
end;

# ╔═╡ 5ca7ac14-1df9-11eb-13f3-e5c86333ef83
begin
	solarSlider = @bind S Slider(1200:2.:1850, default=Sneo);
	md"1200 W/m² $(solarSlider) 1850 W/m²"
end

# ╔═╡ 0bbcdf5a-0dba-11eb-3e81-2b075d4f67ea
begin
	md"""
	*Move the slider below to change the solar insolation:* S = $(S) [W/m²]
	"""
end

# ╔═╡ 5cf4ccdc-0d7e-11eb-2683-fd9a72e763f2
md"""
### (Possible) Homework Exercises:
1. What happens if we instead assume the albedo varies continusouly from -10°C to 10°C? Repeat the above analysis with

$\alpha(T) = \begin{cases}
\alpha_{o} & \mbox{if } T < -10°\text{C}\\
\alpha_{o} + \alpha_{i} \frac{(T+10)}{20} & \mbox{if } -10°\text{C}<T<10°\text{C} \\
\alpha_{o} + \alpha_{i} & \mbox{if } T > 10°\text{C}\\
\end{cases}$

2. What happens if we add CO2 to the neoprotorezoic atmosphere? How much CO2 would we need to add (relevant to a modern background of 280 ppm) to un-freeze the snowball?
"""

# ╔═╡ 1dc709e2-1de8-11eb-28da-fd7469ec50c2
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

# ╔═╡ e9ad66b0-0d6b-11eb-26c0-27413c19dd32
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	Pkg.add("LaTeXStrings")
	using LaTeXStrings
	using Plots
	using PlutoUI
	Model = ingredients("1_energy_balance_model.jl");
end;

# ╔═╡ 262fc3c6-1df2-11eb-332d-c1c9561b3710
function α(T; α0=Model.α, αi=0.5)
	if T >= 0.
		return α0
	elseif T < 0.
		return αi
	end
end

# ╔═╡ 872c8f2a-1df1-11eb-3cfc-3dd568926442
function Model.timestep!(ebm)
	ebm.α = α(ebm.T[end])
	append!(ebm.T, ebm.T[end] + ebm.Δt*Model.tendency(ebm));
	append!(ebm.t, ebm.t[end] + ebm.Δt);
end;

# ╔═╡ 7765f834-0db0-11eb-2c46-ef65f4a1fd1d
begin
	ebm = Model.EBM(Tneo, 0., 1., Model.CO2_const)
	ebm.S = Sneo
	
	ntraj = 20
	Ttraj = repeat([NaN], ntraj)
	Straj = repeat([NaN], ntraj)
end;

# ╔═╡ 477732c4-0daf-11eb-1422-cf0f55cd835b
begin
	S
	restart_ebm!(ebm)
	ebm.S = S
	Model.run!(ebm, 500)

	insert!(Straj, 1, copy(S))
	pop!(Straj)

	insert!(Ttraj, 1, copy(ebm.T[end]))
	pop!(Ttraj)
end;

# ╔═╡ b5849000-0d68-11eb-3beb-c575e8d0ce8e
begin
	Svec = 1200.:1.:1850.
	Svec = vcat(Svec, reverse(Svec))
	Tvec = zeros(size(Svec))

	local T_restart = -100.
	for (i, S) = enumerate(Svec)
		ebm = Model.EBM(T_restart, 0., 1., Model.CO2_const);
		ebm.S = S
		Model.run!(ebm, 500.)
		T_restart = ebm.T[end]
		Tvec[i] = deepcopy(T_restart)
	end
end;

# ╔═╡ 0f222836-0d6c-11eb-2ee8-45545da73cfd
begin
	S
	warming_mask = (1:size(Svec)[1]) .< size(Svec)./2
	p = plot(Svec[warming_mask], Tvec[warming_mask], color="blue",lw=3., alpha=0.5, label="cool branch")
	plot!(Svec[.!warming_mask], Tvec[.!warming_mask], color="red", lw=3., alpha=0.5, label="warm branch")
	plot!(legend=:topleft)
	plot!(xlabel="solar insolation S [W/m²]", ylabel="Global temperature T [°C]")
	plot!([Model.S], [Model.T0], marker=:c, label="present day", color="orange", markersize=8)
	plot!([Sneo], [Tneo], marker=:c, label="neoproterozoic", color="lightblue", markersize=8)
	plot_trajectory!(p, reverse(Straj), reverse(Ttraj), lw=9)
end

# ╔═╡ Cell order:
# ╟─05031b60-1df4-11eb-2b61-956e526b3d4a
# ╟─016c1074-1df4-11eb-2da8-578e25d9456b
# ╟─9c118f9a-1df0-11eb-22dd-b14428994076
# ╟─38346e6a-0d98-11eb-280b-f79787a3c788
# ╠═262fc3c6-1df2-11eb-332d-c1c9561b3710
# ╟─a8dcc0fc-1df8-11eb-21fd-1fdebe5dabfc
# ╠═872c8f2a-1df1-11eb-3cfc-3dd568926442
# ╟─29f90d22-0d7c-11eb-3a04-87ccea95d4f1
# ╟─0bbcdf5a-0dba-11eb-3e81-2b075d4f67ea
# ╟─5ca7ac14-1df9-11eb-13f3-e5c86333ef83
# ╠═0f222836-0d6c-11eb-2ee8-45545da73cfd
# ╠═61fc91ec-1df8-11eb-13c1-098c113b46ec
# ╠═7765f834-0db0-11eb-2c46-ef65f4a1fd1d
# ╠═477732c4-0daf-11eb-1422-cf0f55cd835b
# ╠═b5849000-0d68-11eb-3beb-c575e8d0ce8e
# ╠═9b39df12-1df9-11eb-2eb0-f138980be597
# ╠═90ae18dc-0db8-11eb-0d73-c3b7efaef9b0
# ╟─5cf4ccdc-0d7e-11eb-2683-fd9a72e763f2
# ╠═e9ad66b0-0d6b-11eb-26c0-27413c19dd32
# ╟─1dc709e2-1de8-11eb-28da-fd7469ec50c2
