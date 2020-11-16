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
		"OffsetArrays"
	])
	using Plots
	using PlutoUI
	using Images
	using OffsetArrays
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
	Δx = 0.02
	Δy = 0.02
	Δt = 0.002
	
	κ = 0.01
	
	x = (0. -Δx/2.:Δx:1. +Δx/2.)'
	y = (-1. -Δy/2.:Δy:1. +Δx/2.)
	
	Nx = size(x, 2)
	Ny = size(y, 1)
end;

# ╔═╡ 9036dc6a-204e-11eb-305d-45e760e62bef
begin
	diff_kernel = OffsetArray(zeros(Float64, 3,3), -1:1, -1:1)
	diff_kernel[0, 0] = -4
	diff_kernel[-1, 0] = 1.; diff_kernel[1, 0] = 1.;
	diff_kernel[0, -1] = 1.; diff_kernel[0, 1] = 1.;
	diff_kernel
end

# ╔═╡ 79a0086c-2050-11eb-1974-49d430b5eecd
begin
	function diffuse(T, j, i)
		return κ.*sum(diff_kernel[-1:1,-1:1].*T[j-1:j+1, i-1:i+1])/(2Δx^2)
	end
	diffuse(T) = [diffuse(T, j, i) for j=2:Ny-1, i=2:Nx-1]
end

# ╔═╡ 1cea2b90-205d-11eb-0d06-7df64faf1b53
begin
	adv_kernel = OffsetArray(zeros(Float64, 3,3), -1:1, -1:1)
	adv_kernel[-1, 0] = -1.; adv_kernel[1, 0] = 1.;
	adv_kernel[0, -1] = -1.; adv_kernel[0, 1] = 1.;
	adv_kernel
end

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

# ╔═╡ c4424838-12e2-11eb-25eb-058344b39c8b
begin
	# Initial conditions
	T = -repeat(y, 1, size(x,2));
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
@bind go Clock()

# ╔═╡ 3b0e16a2-12e5-11eb-3130-c763c1c85182


# ╔═╡ 1528ed7e-12e5-11eb-34cf-112d2baa7353
function temperature_heatmap(T)
	p = contourf(x', y, T, color=:bluesreds, levels=-1.25:0.25:1.25, colorbar_title="Temperature [°C]")
	plot!(clims=(-1.25, 1.25))
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

# ╔═╡ 87bfc240-12e3-11eb-03cc-756dc00efa6c
function timestep!(t, T)
	update_ghostcells!(T)
	T[2:end-1, 2:end-1] .+= Δt*(advect(T) .+ diffuse(T))
	t[] += Δt
end;

# ╔═╡ 3b4e4722-12fe-11eb-238d-17aea2c23f58
begin
	CFL_adv = maximum(V)*Δt/Δx
	CFL_diff = κ*Δt/(Δx^2)
	CFL_adv, CFL_diff
end

# ╔═╡ c0e46442-27fb-11eb-2c94-15edbda3f84d
function plot_state()
	X = repeat(xitp(x), size(yitp(y),1), 1)
	Y = repeat(yitp(y), 1, size(xitp(x),2))
	p = temperature_heatmap(T)
	Nq = 15
	quiver!(p, X[(Nq+1)÷2:Nq:end], Y[(Nq+1)÷2:Nq:end], quiver=(U[(Nq+1)÷2:Nq:end]./10., V[(Nq+1)÷2:Nq:end]./10.), color=:black, alpha=0.7)
	plot!(p, xlims=(0., 1.), ylims=(-1.0, 1.0))
	plot!(p, xlabel="longitudinal distance", ylabel="latitudinal distance")
	plot!(p, clabel="Temperature")
	as_png(p)
end

# ╔═╡ bd879bbe-12de-11eb-0d1d-93bba42b6ff9
begin
	go
	nT = 5
	for i = 1:nT
		timestep!(t, T)
	end
	plot_state()
end

# ╔═╡ 3cc1218e-1307-11eb-1907-e7cd68f6af35
heatmap(x', y, ψ̂)

# ╔═╡ d96c7a56-12e4-11eb-123c-d57487bd37df
as_svg(x) = PlutoUI.Show(MIME"image/svg+xml"(), repr(MIME"image/svg+xml"(), x))

# ╔═╡ Cell order:
# ╟─0f8db6f4-2113-11eb-18b4-21a469c67f3a
# ╟─ed741ec6-1f75-11eb-03be-ad6284abaab8
# ╟─ac759b96-2114-11eb-24cb-d50b556f4142
# ╟─3a4a1aea-2118-11eb-30a9-57b87f2ddfae
# ╟─a60e5550-211a-11eb-3cf8-f9bae0a9efd3
# ╠═b1b5625e-211a-11eb-3ee1-3ba9c9cc375a
# ╠═65da5b38-12dc-11eb-3505-bdaf7834afaa
# ╠═9036dc6a-204e-11eb-305d-45e760e62bef
# ╠═fd07ee24-2067-11eb-0ac8-7b3da3993223
# ╠═79a0086c-2050-11eb-1974-49d430b5eecd
# ╠═1cea2b90-205d-11eb-0d06-7df64faf1b53
# ╠═dab0f406-2067-11eb-176d-9dab6819dc98
# ╠═6b3b6030-2066-11eb-3343-e19284638efb
# ╠═16b72cfc-2114-11eb-257d-b7747a99e155
# ╠═b68ca886-2053-11eb-2e39-35c724ed3a3c
# ╠═c4424838-12e2-11eb-25eb-058344b39c8b
# ╠═3b4e4722-12fe-11eb-238d-17aea2c23f58
# ╟─f5ae1756-12e9-11eb-1228-8f03879c154a
# ╟─f9824610-12e7-11eb-3e61-f96c900a0636
# ╠═87bfc240-12e3-11eb-03cc-756dc00efa6c
# ╠═440fe49a-12e5-11eb-1c08-f706f5f33c84
# ╠═bd879bbe-12de-11eb-0d1d-93bba42b6ff9
# ╠═c0e46442-27fb-11eb-2c94-15edbda3f84d
# ╠═3cc1218e-1307-11eb-1907-e7cd68f6af35
# ╠═3b0e16a2-12e5-11eb-3130-c763c1c85182
# ╠═1528ed7e-12e5-11eb-34cf-112d2baa7353
# ╟─bb084ace-12e2-11eb-2dfc-111e90eabfdd
# ╠═627eb1a4-12e2-11eb-30d1-c1ad292d1522
# ╠═e3ee80c0-12dd-11eb-110a-c336bb978c51
# ╠═9c8a7e5a-12dd-11eb-1b99-cd1d52aefa1d
# ╠═d96c7a56-12e4-11eb-123c-d57487bd37df
