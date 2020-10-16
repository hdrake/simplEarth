### A Pluto.jl notebook ###
# v0.12.4

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

# ╔═╡ 232c24de-0e30-11eb-1544-53cb0b45ec9b
begin
	using Pkg
	Pkg.add("Colors")
	Pkg.add("Plots")
	using Colors
	using Plots
	using PlutoUI
	gr()
end

# ╔═╡ beed7b48-0e2b-11eb-2805-edf0be651fc3
md"""
### Numerically solving partial differential equations
#### Learning through one-dimensional advection and diffusion

##### 1. Advection

Consider the differential equation for temperature $T(x,t)$,

$\frac{\partial T}{\partial t} + U \frac{\partial T}{\partial x} = 0$

Since this differential equation derivatives in both time $t$ **and** space $x$, we replace the *ordinary* deriative symbole $d$ with a *partial* derivative symbol $\partial$.

As in our soluation of the "zero-dimenisonal" energy balance model **ODE**, we begin by discretizing time using a *forward finite difference* scheme:

$\frac{T_{n+1}-T_{n}}{\Delta t} + U \frac{\partial T}{\partial x} = 0$
"""

# ╔═╡ c3f99436-0f49-11eb-33b6-e78c150d4ff5
md"""
We now need to discretize the spatial derivative. 

$\frac{T_{n+1,\, i}-T_{n,\, i}}{\Delta t} + U \frac{T_{n,\, i+1} - T_{n,\, i-1}}{2 \Delta x} = 0$

Just as before, we re-order the equation to solve for the temperature of the $i$th grid cell at the next timestep $n+1$,

$T_{n+1,\, i} = T_{n,\, i} + \Delta t \left( U \frac{T_{n,\, i+1} - T_{n,\, i-1}}{2 \Delta x} \right)$

Consider solving this on a grid $i \in [1,\, N_{i}]$. How do we handle the edge case $i=1$, where $T_{n+1,\, 1}$ depends on $T_{n+1,\, 0}$, which is undefined?

For both this point $i=1$ and the other extreme $i=N_{i}$, we need a **boundary condition**. There are a number of ways of doing this, and we will explore some alternative kinds of boundary conditions later, but the simplest is a *periodic boundary condition*, which wraps the grid around by setting $T_{n,\, 0} = T_{n,\, N_{i}}$ and $T_{n,\, 1} = T_{n,\, N_{i}+1}$
"""

# ╔═╡ 0967dc46-0f4d-11eb-13b2-2b2042595b6e
md"### Numerical implementation"

# ╔═╡ 19273c8a-0f4d-11eb-3caf-678f073c6b1a
md"##### Setting up the model parameters

To keep things simple, let us consider our temperature equation as a model of temperature variations in a one-dimensional ocean current of length $L=1$ m and with a speed $U = 1$ m/s.

The choice of the **discretization resolution** $N_{i}$ is up to the modeller, but here we make it $N_{i} = 10$ because we will be able to easily pick out each of the individual grid cells in plots below.
"

# ╔═╡ e700c752-0e2d-11eb-2c61-3959bdba269a
begin
	nx = 10
	Lx = 1.
	Δx = Lx/nx
	Δt = 0.001
	U = 1.
	
	x = Δx/2.:Δx:Lx
end;

# ╔═╡ 5b1b6802-0e2e-11eb-1fa6-1ddd5478cc54
begin
	# Initial conditions
	T = sin.(2π*x);
	t = [0.]
end;

# ╔═╡ 4463cf8c-0ef8-11eb-2c6e-45d982d5687a
function advect!(T)
	T .+= Δt*U*(circshift(T, (1)) .- circshift(T, (-1)))/(2Δx)
end;

# ╔═╡ 67cb0ed2-0efa-11eb-2645-495dba94ea64
function timestep!(t, T)
	advect!(T)
	t .+= Δt
end

# ╔═╡ ba1ea938-0f44-11eb-0fac-cd1afe507d57
@bind go Button("Timestep")

# ╔═╡ 18a3f9dc-0e6d-11eb-0b31-0b5296c2e83b
begin
	⏩ = nothing
	go
	nT = 50
	for i = 1:nT
		timestep!(t, T)
	end
end;

# ╔═╡ f2c7638a-0e34-11eb-2210-9f9b0c3519fe
function temperature_heatmap(T)
	p = plot(xticks=x, yticks=nothing, size=(700,90))
	plot!(p, x, [0.], reshape(T, (size(T)...,1))', st=:heatmap, clims=(-1., 1.))
end;

# ╔═╡ b19bfb58-0f4e-11eb-218f-8d930b7afff6
md"#### 2. Diffusion"

# ╔═╡ d2c101c8-0f4e-11eb-0c4f-0972be924ce1
md"#### 3. Advection-Diffusion"

# ╔═╡ 70ae3138-0f44-11eb-1a8d-d3b5a39c1b42
md"""#### Pluto Book-keeping"""

# ╔═╡ 845e102a-0f44-11eb-3935-2d9a0efecedc
as_svg(x) = PlutoUI.Show(MIME"image/svg+xml"(), repr(MIME"image/svg+xml"(), x))

# ╔═╡ 00cc530a-0e35-11eb-133b-ef8e30ea7b12
begin
	⏩
	p1 = plot(x, T, label="Temperature", ylim=[-1.1, 1.1], xlim=[0., 1.], marker=:c)
	annotate!(p1, [(0.05, 0.9, string("t = ", round(t[1], digits=2)))])
	p2 = temperature_heatmap(T)
	p = plot(p1, p2, layout=grid(2, 1, heights=[0.7 , 0.3]), size=(680,250))
end |> as_svg

# ╔═╡ Cell order:
# ╠═beed7b48-0e2b-11eb-2805-edf0be651fc3
# ╟─c3f99436-0f49-11eb-33b6-e78c150d4ff5
# ╟─0967dc46-0f4d-11eb-13b2-2b2042595b6e
# ╟─19273c8a-0f4d-11eb-3caf-678f073c6b1a
# ╠═e700c752-0e2d-11eb-2c61-3959bdba269a
# ╠═5b1b6802-0e2e-11eb-1fa6-1ddd5478cc54
# ╠═4463cf8c-0ef8-11eb-2c6e-45d982d5687a
# ╠═67cb0ed2-0efa-11eb-2645-495dba94ea64
# ╟─ba1ea938-0f44-11eb-0fac-cd1afe507d57
# ╠═00cc530a-0e35-11eb-133b-ef8e30ea7b12
# ╠═18a3f9dc-0e6d-11eb-0b31-0b5296c2e83b
# ╠═f2c7638a-0e34-11eb-2210-9f9b0c3519fe
# ╟─b19bfb58-0f4e-11eb-218f-8d930b7afff6
# ╟─d2c101c8-0f4e-11eb-0c4f-0972be924ce1
# ╟─70ae3138-0f44-11eb-1a8d-d3b5a39c1b42
# ╠═232c24de-0e30-11eb-1544-53cb0b45ec9b
# ╠═845e102a-0f44-11eb-3935-2d9a0efecedc
