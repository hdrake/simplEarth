### A Pluto.jl notebook ###
# v0.12.10

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

# ╔═╡ 9c8a7e5a-12dd-11eb-1b99-cd1d52aefa1d
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		"Plots",
		"PlutoUI",
		"Images",
		"FileIO",
		"ImageMagick",
		"ImageIO",
		"OffsetArrays",
		"ThreadsX",
		"Strided",
	])
	using Plots
	using PlutoUI
	using Images
	using OffsetArrays
	using ThreadsX
	using Strided
end

# ╔═╡ 0f8db6f4-2113-11eb-18b4-21a469c67f3a
md"""
### Lecture 23: Solving Partial Differential Equations (PDEs) Numerically
**Part II: Heat transport by ocean currents (two-dimensional advection and diffusion)**
"""

# ╔═╡ ed741ec6-1f75-11eb-03be-ad6284abaab8
html"""
<iframe width="700" height="394" src="https://www.youtube-nocookie.com/embed/H4HUJs6LQfI" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ ac759b96-2114-11eb-24cb-d50b556f4142
md"""
### 1) Background: two-dimensional advection-diffusion

##### 1.1) The two-dimensional advection-diffusion equation
Recall from **Lecture 22** that the one-dimensional advection-diffusion equation is written as

$\frac{\partial T(x,t)}{\partial t} = -U \frac{\partial T}{\partial x} + \kappa \frac{\partial^{2} T}{\partial x^{2}},$

where $T(x, t)$ is the temperature, $U$ is a constant advective velocity and $\kappa$ is the diffusivity.

The two-dimensional advection diffusion equation simply adds advection and diffusion operators acting in a second dimensions $y$ (orthogonal to $x$). 

$\frac{\partial T(x,y,t)}{\partial t} = u(x,y) \frac{\partial T}{\partial x} + v(x,y) \frac{\partial T}{\partial y} + \kappa \left( \frac{\partial^{2} T}{\partial x^{2}} + \frac{\partial^{2} T}{\partial y^{2}} \right),$

where $\vec{u}(x,y) = (u, v) = u\,\mathbf{\hat{x}} + v\,\mathbf{\hat{y}}$ is a velocity vector field.

Throughout the rest of the Climate Modelling module, we will consider $x$ to be the *longitundinal* direction (positive from west to east) and $y$ to the be the *latitudinal* direction (positive from south to north).
"""

# ╔═╡ 3a4a1aea-2118-11eb-30a9-57b87f2ddfae
md"""
##### 1.2) Multivariable shorthand notation

Conventionally, the two-dimensional advection-diffusion equation is written more succintly as

$\frac{\partial T(x,y,t)}{\partial t} = - \vec{u} \cdot \nabla T + \kappa \nabla^{2} T,$

using the following shorthand notation.

The **gradient** operator is defined as 

$\nabla \equiv (\frac{\partial}{\partial x}, \frac{\partial }{\partial y})$

such that

$\nabla T = (\frac{\partial T}{\partial x}, \frac{\partial T}{\partial y})$ and 

$\vec{u} \cdot \nabla T = (u, v) \cdot (\frac{\partial T}{\partial x}, \frac{\partial T}{\partial y}) = u \frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y}.$

The **Laplacian** operator $\nabla^{2}$ (sometimes denoted $\Delta$) is defined as 

$\nabla^{2} = \frac{\partial^{2}}{\partial x^{2}} + \frac{\partial^{2}}{\partial y^{2}}$

such that

$\nabla^{2} T = \frac{\partial^{2} T}{\partial x^{2}} + \frac{\partial^{2} T}{\partial y^{2}}.$

The **divergence** operator is defined as $\nabla \cdot [\quad]$, such that

$\nabla \cdot \vec{u} = \left(\frac{\partial}{\partial x}, \frac{\partial}{\partial x} \right) \cdot (u,v) = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}.$

**Note:** Since seawater is largely incompressible, we can approximate ocean currents as a *non-divergent flow*, with $\nabla \cdot \vec{u} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$. Among other implications, this allows us to write:

\begin{align}
\vec{u} \cdot \nabla T&=
u\frac{\partial T(x,y,t)}{\partial x} + v\frac{\partial T(x,y,t)}{\partial y}\newline &=
u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y} + T\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right)\newline &=
\left( u\frac{\partial T}{\partial x} + T\frac{\partial u}{\partial x} \right) +
\left( v\frac{\partial T}{\partial y} + \frac{\partial v}{\partial y} \right)
\newline &=
\frac{\partial (uT)}{\partial x} + \frac{\partial (vT)}{\partial x}\newline &=
\nabla \cdot (\vec{u}T)
\end{align}

using the product rule (separately in both $x$ and $y$).
"""

# ╔═╡ a60e5550-211a-11eb-3cf8-f9bae0a9efd3
md"""
##### 1.3) The flux-form two-dimensional advection-diffusion equation

This lets us finally re-write the two-dimensional advection-diffusion equation as:

$\frac{\partial T}{\partial t} = - \nabla \cdot (\vec{u}T) + \kappa \nabla^{2} T$

which is the form we will use in our numerical algorithm below.
"""

# ╔═╡ b1b5625e-211a-11eb-3ee1-3ba9c9cc375a
md"""
### 2)
"""

# ╔═╡ 65da5b38-12dc-11eb-3505-bdaf7834afaa
begin
	Δx = 0.04
	Δy = 0.04
	
	xs = (0. -Δx/2.:Δx:1. +Δx/2.)'
	ys = (-1. -Δy/2.:Δy:1. +Δx/2.)
	
	Nx = length(xs)
	Ny = length(ys)
end;

# ╔═╡ 490320c0-2818-11eb-1b72-f3c08c502e51
Δt = 0.005

# ╔═╡ 9036dc6a-204e-11eb-305d-45e760e62bef
begin
	diff_kernel = OffsetArray(zeros(Float64, 3,3), -1:1, -1:1)
	diff_kernel[0, 0] = -4
	diff_kernel[-1, 0] = 1.; diff_kernel[1, 0] = 1.;
	diff_kernel[0, -1] = 1.; diff_kernel[0, 1] = 1.;
	diff_kernel
end

# ╔═╡ 1cea2b90-205d-11eb-0d06-7df64faf1b53
begin
	adv_kernel = OffsetArray(zeros(Float64, 3,3), -1:1, -1:1)
	adv_kernel[-1, 0] = -1.; adv_kernel[1, 0] = 1.;
	adv_kernel[0, -1] = -1.; adv_kernel[0, 1] = 1.;
	adv_kernel
end

# ╔═╡ d82eefe0-280e-11eb-1d94-4b2d95630f3c
xs

# ╔═╡ 6b3b6030-2066-11eb-3343-e19284638efb
plot_kernel(A) = heatmap(
	collect(A),
	color=:bluesreds, clims=(-maximum(abs.(A)), maximum(abs.(A))), colorbar=false,
	xticks=false, yticks=false, size=(100, 100), xaxis=false, yaxis=false
)

# ╔═╡ fd07ee24-2067-11eb-0ac8-7b3da3993223
plot_kernel(diff_kernel)

# ╔═╡ dab0f406-2067-11eb-176d-9dab6819dc98
plot_kernel(adv_kernel)

# ╔═╡ b68ca886-2053-11eb-2e39-35c724ed3a3c
function update_ghostcells!(A; option="no-flux")
	Atmp = @view A[:,:]
	if option=="no-flux"
		A[1, :] = Atmp[2, :]; Atmp[end, :] = Atmp[end-1, :]
		A[:, 1] = Atmp[:, 2]; Atmp[:, end] = Atmp[:, end-1]
	end
end

# ╔═╡ 60e62962-280b-11eb-063c-e99186efb59c
function α(T; α0=α0, αi=αi, ΔT=3.)
	if T < -ΔT
		return αi
	elseif -ΔT <= T < ΔT
		return αi + (α0-αi)*(T+ΔT)/(2ΔT)
	elseif T >= ΔT
		return α0
	end
end

# ╔═╡ e0ed2f4a-2836-11eb-20e1-23f8567a435c
plot(-20:20, α, ylim=(0,1))

# ╔═╡ 09c49990-280c-11eb-3ad6-9dd0dd88376d
heatmap(S) |> as_png

# ╔═╡ 96ed4c10-280a-11eb-0eba-139ca63247e3
T0 = 0.0

# ╔═╡ c4424838-12e2-11eb-25eb-058344b39c8b
begin
	# Initial conditions
	T = [
		T0
		for y in ys, x in xs[:]
	]
	t = Ref(0.)
end;

# ╔═╡ f5ae1756-12e9-11eb-1228-8f03879c154a
md"""
### Two-dimensional advection and diffusion
"""

# ╔═╡ f9824610-12e7-11eb-3e61-f96c900a0636
md"""
##### Need boundary conditions still! 
"""

# ╔═╡ 440fe49a-12e5-11eb-1c08-f706f5f33c84
@bind go Clock(.1)

# ╔═╡ 9fe89d82-2833-11eb-2a39-5f94dd51d1ef
@bind temperature_control_event html"""
<script>

const add_button = html`<button>Submit</button>`
const add_field = html`<input type=number value=-5 style='width: 4em;'>`

const set_button = html`<button>Submit</button>`
const set_field = html`<input type=number value=14 style='width: 4em;'>`

const node = html`<div style='border: 1em solid #eeffee; padding: 1em; border-radius: 1.5em;'>
<p>Increase global temperature by ${add_field}°C &nbsp&nbsp&nbsp ${add_button}</p>
<p>Set global temperature to ${set_field}°C &nbsp&nbsp&nbsp ${set_button}</p>

</div>`


add_field.oninput = set_field.oninput = (e) => {
	e.stopPropagation()
}
add_button.onclick = () => {
	node.value = {
		type: "add",
		value: add_field.valueAsNumber
	}
	node.dispatchEvent(new CustomEvent("input"))
}
set_button.onclick = () => {
	node.value = {
		type: "set",
		value: set_field.valueAsNumber
	}
	node.dispatchEvent(new CustomEvent("input"))
}

return node
</script>
"""

# ╔═╡ 4a709370-2818-11eb-1302-0f877986e7b6
κ = 0.015

# ╔═╡ 79a0086c-2050-11eb-1974-49d430b5eecd
begin
	function diffuse(T, j, i)
		return κ.*sum(diff_kernel[-1:1,-1:1].*T[j-1:j+1, i-1:i+1])/(2Δx^2)
	end
	diffuse(T) = [diffuse(T, j, i) for j=2:Ny-1, i=2:Nx-1]
end

# ╔═╡ 9b9e1b30-2810-11eb-2493-717b3949a3c5
diffuse(T)

# ╔═╡ 396ad562-2837-11eb-0a65-55ece65a7da6
α0=0.3

# ╔═╡ 439f7984-2837-11eb-2e08-75d8e90d4b9f
αi=0.5

# ╔═╡ 2ae5d9ee-280b-11eb-1bfc-c79d2742eee8
A = 11

# ╔═╡ 30539f30-280b-11eb-1219-9b1d7d8c8cfb
B = -0.7

# ╔═╡ 13f114d0-280b-11eb-2c89-67cc659372ce
function outgoing_thermal_radiation(T)
	A .- B .* (T .- T0)
end

# ╔═╡ 5f6ca7c2-2816-11eb-2bdb-6502b84cc9c8
plot(
	-10:40, outgoing_thermal_radiation(-10:40),
	xlabel="Temperature",
	ylabel="Outgoing radiation"
)

# ╔═╡ 06367924-2839-11eb-1de9-51833725a659
S_peak = 35.0

# ╔═╡ ec3798f0-280b-11eb-3e26-9d40d35a6920
Ss = [
		S_peak .* (cos((y+1) * π/4) + .3)
		for y in ys, x in xs[:]
	]

# ╔═╡ 443fb830-280b-11eb-017e-0d8a56cf7729
function absorbed_solar_radiation(T)
	absorption = 1.0 .- α.(T)
	
	absorption .* Ss
end

# ╔═╡ b2c066b0-2815-11eb-0373-453a84ff3d3b
mean(x) = sum(x) / length(x)

# ╔═╡ bff76d60-2815-11eb-2dc3-8f92794a5056
go; mean(T)

# ╔═╡ b0a94630-2833-11eb-1f37-63e7f5beaf10
begin
	if !ismissing(temperature_control_event)
		e = temperature_control_event
		if e["type"] == "add"
			T .+= e["value"]
		elseif e["type"] == "set"
			T .= e["value"]
		end
		
	end
	temperature_control_event_handled = true
	Text("Temperature control logic")
end

# ╔═╡ 3b0e16a2-12e5-11eb-3130-c763c1c85182


# ╔═╡ 1528ed7e-12e5-11eb-34cf-112d2baa7353
function temperature_heatmap(T)
	levels = -10:1.0:40
	
	p = contourf(xs', ys, T, 
		levels=levels, 
		color=:bluesreds, colorbar_title="Temperature [°C]",
		lw=0,
		clims=extrema(levels)
	)
	
	contour!(p, xs', ys, T, 
		levels=[0], 
		color=:white,
		lw=3
	)
end

# ╔═╡ bb084ace-12e2-11eb-2dfc-111e90eabfdd
md"##### Setting up the velocity field"

# ╔═╡ e3ee80c0-12dd-11eb-110a-c336bb978c51
begin
	∂x(ϕ) = (ϕ[:,2:end] - ϕ[:,1:end-1])/Δx
	∂y(ϕ) = (ϕ[2:end,:] - ϕ[1:end-1,:])/Δy
	
	xpad(ϕ) = hcat(zeros(size(ϕ,1)), ϕ, zeros(size(ϕ,1)))
	ypad(ϕ) = vcat(zeros(size(ϕ,2))', ϕ, zeros(size(ϕ,2))')
	
	xitp(ϕ) = 0.5*(ϕ[:,2:end]+ϕ[:,1:end-1])
	yitp(ϕ) = 0.5*(ϕ[2:end,:]+ϕ[1:end-1,:])
	
	function diagnose_velocities(ψ)
		u = ∂y(ψ)
		v = -∂x(ψ)
		return u,v
	end
end

# ╔═╡ 627eb1a4-12e2-11eb-30d1-c1ad292d1522
begin
	ϵ = 0.05
	xψ = (-Δx:Δx:1. +Δx)'
	yψ = (-1-Δy:Δy:1. +Δy)
	
	# See page 595 of Vallis Edt.2
	ψ̂(x,y) = π*sin.(π*y) * (
		1 .- x - exp.(-x/(2*ϵ)) .* (
			cos.(√3*x/(2*ϵ)) .+
			(1. /√3)*sin.(√3*x/(2*ϵ))
			)
		.+ ϵ*exp.((x .- 1.)/ϵ)
	)
	
	u,v = diagnose_velocities(ψ̂(xψ, yψ))
	U = xitp(u) ./10.
	V = yitp(v) ./10.
	U[1,:] .= 0.; V[1,:] .= 0.;
	U[end,:] .= 0.; V[end,:] .= 0.;
	U[:,1] .= 0.; V[:,1] .= 0.;
	U[:,end] .= 0.; V[:,end] .= 0.;
end;

# ╔═╡ 16b72cfc-2114-11eb-257d-b7747a99e155
begin
	function advect(T, j, i)
		return .-(
			sum(adv_kernel[0, -1:1].*(U[j, i-1:i+1].*T[j, i-1:i+1]))/(2Δx) .+
			sum(adv_kernel[-1:1, 0].*(V[j-1:j+1, i].*T[j-1:j+1, i]))/(2Δy)
		)
	end
	advect(T) = [advect(T, j, i) for j=2:Ny-1, i=2:Nx-1]
end

# ╔═╡ 918ae9c0-2810-11eb-0620-956ef6dac50c
advect(T)

# ╔═╡ 42ea047e-2811-11eb-0e60-510310f813d8
function timestep2!(t, T)
	update_ghostcells!(T)
	T[2:end-1, 2:end-1] .+= Δt.*(
		advect(T) .+ diffuse(T) .+ 
		(absorbed_solar_radiation(T)[2:end-1, 2:end-1]) .- 
		(outgoing_thermal_radiation(T)[2:end-1, 2:end-1])
	)
	t[] += Δt
end;

# ╔═╡ ad1aec70-2811-11eb-0473-bbb02625902a
for i in 1:10
	timestep2!(t, T)
end

# ╔═╡ 87bfc240-12e3-11eb-03cc-756dc00efa6c
function timestep!(t, T)
	update_ghostcells!(T)
	T[2:end-1, 2:end-1] .+= Δt*(
		advect(T) .+ diffuse(T) .+ 
		(@view absorbed_solar_radiation(T)[2:end-1, 2:end-1]) .- 
		(@view outgoing_thermal_radiation(T)[2:end-1, 2:end-1])
	)
	t[] += Δt
end;

# ╔═╡ 1a880910-2811-11eb-15db-158d2c8eb1a4
for i in 1:10
	timestep!(t, T)
end

# ╔═╡ 3b4e4722-12fe-11eb-238d-17aea2c23f58
begin
	CFL_adv = maximum(V)*Δt/Δx
	CFL_diff = κ*Δt/(Δx^2)
	CFL_adv, CFL_diff
end

# ╔═╡ c0e46442-27fb-11eb-2c94-15edbda3f84d
function plot_state()
	X = repeat(xitp(xs), size(yitp(ys),1), 1)
	Y = repeat(yitp(ys), 1, size(xitp(xs),2))
	p = temperature_heatmap(T)
	Nq = 4
	quiver!(p, X[(Nq+1)÷2:Nq:end], Y[(Nq+1)÷2:Nq:end], quiver=(U[(Nq+1)÷2:Nq:end]./10., V[(Nq+1)÷2:Nq:end]./10.), color=:black, alpha=0.7)
	plot!(p, xlims=(0., 1.), ylims=(-1.0, 1.0))
	plot!(p, xlabel="longitudinal distance", ylabel="latitudinal distance")
	plot!(p, clabel="Temperature")
	as_png(p)
end

# ╔═╡ bd879bbe-12de-11eb-0d1d-93bba42b6ff9
begin
	go
	temperature_control_event_handled
	nT = 50
	for i = 1:nT
		timestep!(t, T)
	end
	plot_state()
end

# ╔═╡ 3cc1218e-1307-11eb-1907-e7cd68f6af35
heatmap(xs', ys, ψ̂)

# ╔═╡ Cell order:
# ╟─0f8db6f4-2113-11eb-18b4-21a469c67f3a
# ╟─ed741ec6-1f75-11eb-03be-ad6284abaab8
# ╟─ac759b96-2114-11eb-24cb-d50b556f4142
# ╟─3a4a1aea-2118-11eb-30a9-57b87f2ddfae
# ╟─a60e5550-211a-11eb-3cf8-f9bae0a9efd3
# ╠═b1b5625e-211a-11eb-3ee1-3ba9c9cc375a
# ╠═65da5b38-12dc-11eb-3505-bdaf7834afaa
# ╠═490320c0-2818-11eb-1b72-f3c08c502e51
# ╠═9036dc6a-204e-11eb-305d-45e760e62bef
# ╠═fd07ee24-2067-11eb-0ac8-7b3da3993223
# ╠═79a0086c-2050-11eb-1974-49d430b5eecd
# ╠═1cea2b90-205d-11eb-0d06-7df64faf1b53
# ╠═dab0f406-2067-11eb-176d-9dab6819dc98
# ╠═d82eefe0-280e-11eb-1d94-4b2d95630f3c
# ╠═6b3b6030-2066-11eb-3343-e19284638efb
# ╠═16b72cfc-2114-11eb-257d-b7747a99e155
# ╠═b68ca886-2053-11eb-2e39-35c724ed3a3c
# ╠═918ae9c0-2810-11eb-0620-956ef6dac50c
# ╠═9b9e1b30-2810-11eb-2493-717b3949a3c5
# ╠═13f114d0-280b-11eb-2c89-67cc659372ce
# ╠═e0ed2f4a-2836-11eb-20e1-23f8567a435c
# ╠═60e62962-280b-11eb-063c-e99186efb59c
# ╠═ec3798f0-280b-11eb-3e26-9d40d35a6920
# ╠═09c49990-280c-11eb-3ad6-9dd0dd88376d
# ╠═443fb830-280b-11eb-017e-0d8a56cf7729
# ╠═c4424838-12e2-11eb-25eb-058344b39c8b
# ╠═96ed4c10-280a-11eb-0eba-139ca63247e3
# ╠═3b4e4722-12fe-11eb-238d-17aea2c23f58
# ╠═1a880910-2811-11eb-15db-158d2c8eb1a4
# ╠═ad1aec70-2811-11eb-0473-bbb02625902a
# ╠═42ea047e-2811-11eb-0e60-510310f813d8
# ╟─f5ae1756-12e9-11eb-1228-8f03879c154a
# ╟─f9824610-12e7-11eb-3e61-f96c900a0636
# ╠═5f6ca7c2-2816-11eb-2bdb-6502b84cc9c8
# ╠═87bfc240-12e3-11eb-03cc-756dc00efa6c
# ╠═440fe49a-12e5-11eb-1c08-f706f5f33c84
# ╟─9fe89d82-2833-11eb-2a39-5f94dd51d1ef
# ╟─bd879bbe-12de-11eb-0d1d-93bba42b6ff9
# ╠═bff76d60-2815-11eb-2dc3-8f92794a5056
# ╠═4a709370-2818-11eb-1302-0f877986e7b6
# ╠═396ad562-2837-11eb-0a65-55ece65a7da6
# ╠═439f7984-2837-11eb-2e08-75d8e90d4b9f
# ╠═2ae5d9ee-280b-11eb-1bfc-c79d2742eee8
# ╠═30539f30-280b-11eb-1219-9b1d7d8c8cfb
# ╠═06367924-2839-11eb-1de9-51833725a659
# ╠═b2c066b0-2815-11eb-0373-453a84ff3d3b
# ╠═c0e46442-27fb-11eb-2c94-15edbda3f84d
# ╠═b0a94630-2833-11eb-1f37-63e7f5beaf10
# ╠═3cc1218e-1307-11eb-1907-e7cd68f6af35
# ╠═3b0e16a2-12e5-11eb-3130-c763c1c85182
# ╠═1528ed7e-12e5-11eb-34cf-112d2baa7353
# ╟─bb084ace-12e2-11eb-2dfc-111e90eabfdd
# ╠═627eb1a4-12e2-11eb-30d1-c1ad292d1522
# ╠═e3ee80c0-12dd-11eb-110a-c336bb978c51
# ╠═9c8a7e5a-12dd-11eb-1b99-cd1d52aefa1d
