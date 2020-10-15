### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 232c24de-0e30-11eb-1544-53cb0b45ec9b
begin
	using Pkg
	Pkg.add("Colors")
	Pkg.add("Plots")
	using Colors
	using Plots
	gr()
end

# ╔═╡ beed7b48-0e2b-11eb-2805-edf0be651fc3
md"""
### Numerically solving partial differential equations
##### Advection and Diffusion

Consider the differential equation for temperature $T(x,t)$,

$\frac{\partial T}{\partial t} + U \frac{\partial T}{\partial x} = 0$

Since this differential equation derivatives in both time $t$ **and** space $x$, we replace the *ordinary* deriative symbole $d$ with a *partial* derivative symbol $\partial$.

As in our soluation of the "zero-dimenisonal" energy balance model **ODE**, we begin by discretizing time using a *forward finite difference* scheme:

$\frac{T_{n+1}-T_{n}}{\Delta t} + U \frac{\partial T}{\partial x} = 0$

We now need to discretize the spatial derivative. This is a bit more ambiguous than it was for the time dimension $t$, where it is obvious that we want to advance *forward* in time. For the moment, let us just do a *right-side forward difference* (towards increasing $x$).

"""

# ╔═╡ e700c752-0e2d-11eb-2c61-3959bdba269a
begin
	nx = 10
	Lx = 1.
	Δx = Lx/nx
	Δt = 0.001
	U = 1.
	
	x = Δx/2.:Δx:Lx
end

# ╔═╡ 254fabac-0e2e-11eb-384d-1b4ac207ccab
initialize_T(x) = exp.(-((x.-Lx/2.)/0.3).^2);

# ╔═╡ 5b1b6802-0e2e-11eb-1fa6-1ddd5478cc54
T = initialize_T(x);

# ╔═╡ 4463cf8c-0ef8-11eb-2c6e-45d982d5687a
function advect!(T)
	T .+= Δt*U*(circshift(T, (1)) .- circshift(T, (-1)))/(2Δx)
end;

# ╔═╡ 18a3f9dc-0e6d-11eb-0b31-0b5296c2e83b
begin
	nT = 30
	for i = 1:nT
		advect!(T)
	end
	⏩ = nothing # dummy variable that triggers plot to udpate
end

# ╔═╡ f2c7638a-0e34-11eb-2210-9f9b0c3519fe
function plot_temperature(T)
	p = plot(xticks=x, yticks=nothing, size=(700,90))
	plot!(p, x, [0.], reshape(T, (size(T)...,1))', st=:heatmap, clims=(0., 1.))
end

# ╔═╡ 00cc530a-0e35-11eb-133b-ef8e30ea7b12
begin
	⏩
	p1 = plot(x, T, label="Temperature", ylim=[0., 1.], xlim=[0., 1.], marker=:c)
	p2 = plot_temperature(T)
	plot(p1, p2, layout=grid(2, 1, heights=[0.7 , 0.3]), size=(680,250))
end

# ╔═╡ Cell order:
# ╠═beed7b48-0e2b-11eb-2805-edf0be651fc3
# ╠═e700c752-0e2d-11eb-2c61-3959bdba269a
# ╠═254fabac-0e2e-11eb-384d-1b4ac207ccab
# ╠═5b1b6802-0e2e-11eb-1fa6-1ddd5478cc54
# ╠═4463cf8c-0ef8-11eb-2c6e-45d982d5687a
# ╟─00cc530a-0e35-11eb-133b-ef8e30ea7b12
# ╠═18a3f9dc-0e6d-11eb-0b31-0b5296c2e83b
# ╠═f2c7638a-0e34-11eb-2210-9f9b0c3519fe
# ╠═232c24de-0e30-11eb-1544-53cb0b45ec9b
