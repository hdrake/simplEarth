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
<iframe width="700" height="394" src="https://www.youtube-nocookie.com/embed/6_GQuVopmUM?start=15" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ ac759b96-2114-11eb-24cb-d50b556f4142
md"""
### 1) Background: two-dimensional advection-diffusion

##### 1.1) The two-dimensional advection-diffusion equation
Recall from **Lecture 22** that the one-dimensional advection-diffusion equation is written as

$\frac{\partial T(x,t)}{\partial t} = -U \frac{\partial T}{\partial x} + \kappa \frac{\partial^{2} T}{\partial x^{2}} - (A - BT)/C + S(1-α)/4/C,$

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

# ╔═╡ cd2ee4ca-2a06-11eb-0e61-e9a2ecf72bd6
struct Grid
	N::Int64
	L::Float64
	
	Δx::Float64
	Δy::Float64
	
	x::Array{Float64, 2}
	y::Array{Float64, 2}
	
	Nx::Int64
	Ny::Int64
	
	function Grid(N, L)
		Δx = L/N # [m]
		Δy = L/N # [m]
		
		x = 0. -Δx/2.:Δx:L +Δx/2.
		x = reshape(x, (1, size(x,1)))
		y = -L -Δy/2.:Δy:L +Δy/2.
		y = reshape(y, (size(y,1), 1))

		Nx, Ny = size(x, 2), size(y, 1)
		
		return new(N, L, Δx, Δy, x, y, Nx, Ny)
	end
end

# ╔═╡ 2a93145e-2a09-11eb-323b-01817062aa89
struct Parameters
	κ::Float64
end

# ╔═╡ 9c3643ea-2a09-11eb-2d38-1589b871ae4a
begin
	import Base.zeros
	zeros(G::Grid) = zeros(G.Ny, G.Nx)
end

# ╔═╡ 32663184-2a81-11eb-0dd1-dd1e10ed9ec6
abstract type ClimateModel
end

# ╔═╡ d3796644-2a05-11eb-11b8-87b6e8c311f9
begin
	struct OceanModel <: ClimateModel
		G::Grid
		P::Parameters

		u::Array{Float64, 2}
		v::Array{Float64, 2}
	end

	OceanModel(G, P) = OceanModel(G, P, zeros(G), zeros(G))
	OceanModel(G) = OceanModel(G, Parameters(1.e4), zeros(G), zeros(G))
end;

# ╔═╡ f92086c4-2a74-11eb-3c72-a1096667183b
begin
	struct ClimateModelSimulation
		C::ClimateModel
		
		T::Array{Float64, 2}
		Δt::Float64
		iter::Base.RefValue{Int64}
	end
	ClimateModelSimulation(C, T, Δt) = ClimateModelSimulation(C, T, Δt, Ref(0))
end

# ╔═╡ 3d12c114-2a0a-11eb-131e-d1a39b4f440b
function InitBox(G; value=1., nx=2, ny=2)
	T = zeros(G)
	T[G.Ny÷2-ny:G.Ny÷2+ny, G.Nx÷2-nx:G.Nx÷2+nx] .= value
	return T
end

# ╔═╡ 2b3d2062-2a73-11eb-368e-f563b2ad7aba
md"""
Let's check that the total heat of the ocean is conserved: 
"""

# ╔═╡ d042d25a-2a62-11eb-33fe-65494bb2fad5
begin
	quiverBox = @bind show_quiver CheckBox(default=false)
	anomalyBox = @bind show_anomaly CheckBox(default=false)
	md"""
	*Click to show the velocity field* $(quiverBox) *or to show temperature **anomalies** instead of absolute values* $(anomalyBox)
	"""
end

# ╔═╡ 933d42fa-2a67-11eb-07de-61cab7567d7d
begin
	κ_ex_Slider = @bind κ_ex Slider(0.:1e5, default=1e4, show_value=true)
	md"""
	*Vary the diffusivity κ* $(κ_ex_Slider)
	"""
end

# ╔═╡ c9ea0f72-2a67-11eb-20ba-376ca9c8014f
@bind go_ex Clock(0.3)

# ╔═╡ 9036dc6a-204e-11eb-305d-45e760e62bef
begin
	xdiff_kernel = OffsetArray(reshape([1., -2., 1.], 1, 3),  0:0, -1:1)
	ydiff_kernel = OffsetArray(reshape([1., -2., 1.], 3, 1),  -1:1, 0:0)
end;

# ╔═╡ 79a0086c-2050-11eb-1974-49d430b5eecd
begin
	function diffuse(T, κ, Δy, Δx, j, i)
		return κ.*(
			sum(xdiff_kernel[0, -1:1].*T[j, i-1:i+1])/(Δx^2) +
			sum(ydiff_kernel[-1:1, 0].*T[j-1:j+1, i])/(Δy^2)
		)
	end
	diffuse(T, κ, Δy, Δx) = [
		diffuse(T, κ, Δy, Δx, j, i) for j=2:size(T, 1)-1, i=2:size(T, 2)-1
	]
	diffuse(T, O::OceanModel) = diffuse(T, O.P.κ, O.G.Δy, O.G.Δx)
end

# ╔═╡ 1cea2b90-205d-11eb-0d06-7df64faf1b53
begin
	adv_kernel = OffsetArray(zeros(Float64, 3,3), -1:1, -1:1)
	adv_kernel[-1, 0] = -1.; adv_kernel[1, 0] = 1.;
	adv_kernel[0, -1] = -1.; adv_kernel[0, 1] = 1.;
	adv_kernel
end

# ╔═╡ 16b72cfc-2114-11eb-257d-b7747a99e155
begin
	function advect(T, u, v, Δy, Δx, j, i)
		return .-(
			sum(adv_kernel[0, -1:1].*(u[j, i-1:i+1].*T[j, i-1:i+1]))/(2Δx) .+
			sum(adv_kernel[-1:1, 0].*(v[j-1:j+1, i].*T[j-1:j+1, i]))/(2Δy)
		)
	end
	advect(T, u, v, Δy, Δx) = [
		advect(T, u, v, Δy, Δx, j, i)
		for j=2:size(T, 1)-1, i=2:size(T, 2)-1
	]
	
	advect(T, O::OceanModel) = advect(T, O.u, O.v, O.G.Δy, O.G.Δx)
end

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
	linearT(G) = 0.5*(1. .+[ -(y/G.L) for y in G.y[:, 1], x in G.x[1, :] ])
	t = Ref(0.)
end;

# ╔═╡ d2c27a5c-2a77-11eb-0772-3b8b0d4cf682
md"""
$\mathcal{O}(\text{timestep})*\mathcal{O}(\text{diffusion / advection}) = N * N^2$
"""

# ╔═╡ a138d27e-2a77-11eb-2d5a-09a717f94e7d
md"$\Delta t < \frac{\Delta x}{max(v)}$"

# ╔═╡ f5ae1756-12e9-11eb-1228-8f03879c154a
md"""
### Two-dimensional advection and diffusion
"""

# ╔═╡ f9824610-12e7-11eb-3e61-f96c900a0636
md"""
##### Need boundary conditions still! 
"""

# ╔═╡ 87bfc240-12e3-11eb-03cc-756dc00efa6c
function timestep!(S::ClimateModelSimulation)
	update_ghostcells!(S.T)
	S.T[2:end-1, 2:end-1] .+= S.Δt*(advect(S.T, S.C) .+ diffuse(S.T, S.C))
	S.iter[] += 1
end;

# ╔═╡ 16905a6a-2a78-11eb-19ea-81adddc21088
begin
	
	Nvec = 1:25
	tvec = map(Nvec) do Npower
		G = Grid(8*Npower, 6.e6);
		P = Parameters(κ_ex);

		#u, v = DoubleGyre(G)
		#u, v = PointVortex(G, Ω=0.5)
		u, v = zeros(G), zeros(G)
		
		O = OceanModel(G, P, u, v)
		
		IC = InitBox(G)
		#IC = InitBox(G, nx=G.Nx÷2-1)
		#IC = linearT(G)

		Δt = 6*60*60
		S = ClimateModelSimulation(O, copy(IC), Δt)

		return @elapsed timestep!(S)
	end
end

# ╔═╡ bb084ace-12e2-11eb-2dfc-111e90eabfdd
md"""##### Computing a quasi-realistic ocean velocity field $\vec{u} = (u, v)$
Our velocity field is given by an analytical solution to the classic wind-driven gyre
problem, which is given by solving the fourth-order partial differential equation:

$- \epsilon_{M} \hat{\nabla}^{4} \hat{\Psi} + \frac{\partial \hat{\Psi} }{ \partial \hat{x}} = \nabla \times \hat{\tau} \mathbf{z},$

where the hats denote that all of the variables have been non-dimensionalized and all of their constant coefficients have been bundles into the single parameter $\epsilon_{M} \equiv \dfrac{\nu}{\beta L^3}$.

The solution makes use of an advanced *asymptotic method* (valid in the limit that $\epsilon \ll 1$) known as *boundary layer analysis* (see MIT course 18.305 to learn more). 
"""



# ╔═╡ df706ebc-2a63-11eb-0b09-fd9f151cb5a8
function impose_no_flux!(u, v)
	u[1,:] .= 0.; v[1,:] .= 0.;
	u[end,:] .= 0.; v[end,:] .= 0.;
	u[:,1] .= 0.; v[:,1] .= 0.;
	u[:,end] .= 0.; v[:,end] .= 0.;
end

# ╔═╡ e3ee80c0-12dd-11eb-110a-c336bb978c51
begin
	∂x(ϕ, Δx) = (ϕ[:,2:end] - ϕ[:,1:end-1])/Δx
	∂y(ϕ, Δy) = (ϕ[2:end,:] - ϕ[1:end-1,:])/Δy
	
	xpad(ϕ) = hcat(zeros(size(ϕ,1)), ϕ, zeros(size(ϕ,1)))
	ypad(ϕ) = vcat(zeros(size(ϕ,2))', ϕ, zeros(size(ϕ,2))')
	
	xitp(ϕ) = 0.5*(ϕ[:,2:end]+ϕ[:,1:end-1])
	yitp(ϕ) = 0.5*(ϕ[2:end,:]+ϕ[1:end-1,:])
	
	function diagnose_velocities(ψ, G)
		u = xitp(∂y(ψ, G.Δy/G.L))
		v = yitp(-∂x(ψ, G.Δx/G.L))
		return u,v
	end
end

# ╔═╡ e2e4cfac-2a63-11eb-1b7f-9d8d5d304b43
function PointVortex(G; Ω=1., a=0.2, x0=0.5, y0=0.)
	x = reshape(0. -G.Δx/(G.L):G.Δx/G.L:1. +G.Δx/(G.L), (1, G.Nx+1))
	y = reshape(-1. -G.Δy/(G.L):G.Δy/G.L:1. +G.Δy/(G.L), (G.Ny+1, 1))
	
	function ψ̂(x,y)
		r = sqrt.((y .-y0).^2 .+ (x .-x0).^2)
		
		stream = -Ω/4*r.^2
		stream[r .> a] = -Ω*a^2/4*(1. .+ 2*log.(r[r .> a]/a))
		
		return stream
	end
		
	u, v = diagnose_velocities(ψ̂(x, y), G)
	impose_no_flux!(u, v)
	
	return u,v
end

# ╔═╡ 863a6330-2a08-11eb-3992-c3db439fb624
begin
	G = Grid(20, 6.e6);
	P = Parameters(κ_ex);
	
	#u, v = DoubleGyre(G)
	u, v = PointVortex(G, Ω=0.5)
	#u, v = zeros(G), zeros(G)

	#IC = InitBox(G)
	IC = InitBox(G, nx=G.Nx÷2-1)
	#IC = linearT(G)
	
	C = OceanModel(G, P, u, v)
	Δt = 6*60*60
	
	S = ClimateModelSimulation(C, copy(IC), Δt)
end;

# ╔═╡ 3b4e4722-12fe-11eb-238d-17aea2c23f58
begin
	CFL_adv(S::ClimateModelSimulation) = maximum(S.C.v)*S.Δt/S.C.G.Δx
	CFL_diff(S::ClimateModelSimulation) = S.C.P.κ*S.Δt/(S.C.G.Δx^2)
	CFL_adv(S), CFL_diff(S)
end

# ╔═╡ c0e46442-27fb-11eb-2c94-15edbda3f84d
function plot_state(C; clims=(-1.1,1.1), show_quiver=true, IC=nothing)
	O = C.C
	X = repeat(O.G.x, O.G.Ny, 1)
	Y = repeat(O.G.y, 1, O.G.Nx)
	if show_anomaly
		arrow_col = :black
		maxdiff = maximum(abs.(O.T .- IC))
		p = contourf(O.G.x[:], O.G.y[:], C.T .- IC, clims=(-maxdiff, maxdiff),
			color=:balance, colorbar_title="Temperature [°C]", linewidth=0.,
			size=(400,530)
		)
	else
		arrow_col = :lightblue
		p = contourf(O.G.x[:], O.G.y[:], C.T,
			color=:thermal, levels=clims[1]:(clims[2]-clims[1])/21.:clims[2],
			colorbar_title="Temperature [°C]", clims=clims,
			linewidth=0., size=(400,530)
		)
	end
	annotate!(500e3, 5500e3,
		text(
			string("t = ", Int64(round(S.iter[]*S.Δt/(60*60*24))), " days"),
			color=:white, :left, 10
		)
	)
	
	Nq = O.G.N÷5
	if show_quiver
		quiver!(p,
			X[(Nq+1)÷2:Nq:end], Y[(Nq+1)÷2:Nq:end],
			quiver=O.G.L*4 .*(O.u[(Nq+1)÷2:Nq:end], O.v[(Nq+1)÷2:Nq:end]),
			color=arrow_col, alpha=0.7
		)
	end
	plot!(p,
		yticks=( (-O.G.L:1000e3:O.G.L), 1e-3*(-O.G.L:1000e3:O.G.L) ),
		xticks=( (0:1000e3:O.G.L), 1e-3*(0:1000e3:O.G.L) ),
		xlims=(0., O.G.L), ylims=(-O.G.L, O.G.L),
	)
	plot!(p, xlabel="longitudinal distance [km]", ylabel="latitudinal distance [km]")
	plot!(p, clabel="Temperature")
	as_png(p)
end

# ╔═╡ ecaab27e-2a16-11eb-0e99-87c91e659cf3
function DoubleGyre(G; β=2e-11, τ₀=0.1, ρ₀=1.e3, ν=1.e5, κ=1.e5, H=1000.)
	ϵM = ν/(β*G.L^3)
	ϵ = ϵM^(1/3.)
	x = reshape(0. -G.Δx/(G.L):G.Δx/G.L:1. +G.Δx/(G.L), (1, G.Nx+1))
	y = reshape(-1. -G.Δy/(G.L):G.Δy/G.L:1. +G.Δy/(G.L), (G.Ny+1, 1))
	
	ψ̂(x,y) = π*sin.(π*y) * (
		1 .- x - exp.(-x/(2*ϵ)) .* (
			cos.(√3*x/(2*ϵ)) .+
			(1. /√3)*sin.(√3*x/(2*ϵ))
			)
		.+ ϵ*exp.((x .- 1.)/ϵ)
	)
		
	u, v = (τ₀/ρ₀)/(β*G.L*H) .* diagnose_velocities(ψ̂(x, y), G)
	impose_no_flux!(u, v)
	
	return u, v
end

# ╔═╡ d96c7a56-12e4-11eb-123c-d57487bd37df
as_svg(x) = PlutoUI.Show(MIME"image/svg+xml"(), repr(MIME"image/svg+xml"(), x))

# ╔═╡ 83c5dbb2-2a0a-11eb-0b1d-d120efa14de5
begin
 	go_ex	
	for i in 1:150
		timestep!(S)
	end
	plot_state(S, clims=(0., 1), show_quiver=show_quiver, IC=IC)
end |> as_svg

# ╔═╡ 794c2148-2a78-11eb-2756-5bd28b7726fa
begin
	plot(8*Nvec, tvec, xlabel="Number of Grid Cells (in x-direction)", ylabel="elapsed time per timestep [s]")
end |> as_svg

# ╔═╡ 6b3b6030-2066-11eb-3343-e19284638efb
plot_kernel(A) = heatmap(
	collect(A),
	color=:bluesreds, clims=(-maximum(abs.(A)), maximum(abs.(A))), colorbar=false,
	xticks=false, yticks=false, size=(30+30*size(A, 2), 30+30*size(A, 1)), xaxis=false, yaxis=false
)

# ╔═╡ fd07ee24-2067-11eb-0ac8-7b3da3993223
plot_kernel(xdiff_kernel[0:0, :]),
plot_kernel(ydiff_kernel[:, 0:0])

# ╔═╡ dab0f406-2067-11eb-176d-9dab6819dc98
plot_kernel(adv_kernel),
plot_kernel(adv_kernel[0:0, :]),
plot_kernel(adv_kernel[:, 0:0])

# ╔═╡ Cell order:
# ╟─0f8db6f4-2113-11eb-18b4-21a469c67f3a
# ╟─ed741ec6-1f75-11eb-03be-ad6284abaab8
# ╟─ac759b96-2114-11eb-24cb-d50b556f4142
# ╟─3a4a1aea-2118-11eb-30a9-57b87f2ddfae
# ╟─a60e5550-211a-11eb-3cf8-f9bae0a9efd3
# ╠═b1b5625e-211a-11eb-3ee1-3ba9c9cc375a
# ╠═cd2ee4ca-2a06-11eb-0e61-e9a2ecf72bd6
# ╠═2a93145e-2a09-11eb-323b-01817062aa89
# ╠═9c3643ea-2a09-11eb-2d38-1589b871ae4a
# ╠═32663184-2a81-11eb-0dd1-dd1e10ed9ec6
# ╠═d3796644-2a05-11eb-11b8-87b6e8c311f9
# ╠═f92086c4-2a74-11eb-3c72-a1096667183b
# ╠═863a6330-2a08-11eb-3992-c3db439fb624
# ╟─3d12c114-2a0a-11eb-131e-d1a39b4f440b
# ╟─2b3d2062-2a73-11eb-368e-f563b2ad7aba
# ╟─d042d25a-2a62-11eb-33fe-65494bb2fad5
# ╠═933d42fa-2a67-11eb-07de-61cab7567d7d
# ╟─c9ea0f72-2a67-11eb-20ba-376ca9c8014f
# ╠═83c5dbb2-2a0a-11eb-0b1d-d120efa14de5
# ╠═9036dc6a-204e-11eb-305d-45e760e62bef
# ╟─fd07ee24-2067-11eb-0ac8-7b3da3993223
# ╠═79a0086c-2050-11eb-1974-49d430b5eecd
# ╠═1cea2b90-205d-11eb-0d06-7df64faf1b53
# ╟─dab0f406-2067-11eb-176d-9dab6819dc98
# ╠═16b72cfc-2114-11eb-257d-b7747a99e155
# ╠═b68ca886-2053-11eb-2e39-35c724ed3a3c
# ╠═c4424838-12e2-11eb-25eb-058344b39c8b
# ╠═16905a6a-2a78-11eb-19ea-81adddc21088
# ╠═794c2148-2a78-11eb-2756-5bd28b7726fa
# ╠═d2c27a5c-2a77-11eb-0772-3b8b0d4cf682
# ╠═a138d27e-2a77-11eb-2d5a-09a717f94e7d
# ╠═3b4e4722-12fe-11eb-238d-17aea2c23f58
# ╟─f5ae1756-12e9-11eb-1228-8f03879c154a
# ╟─f9824610-12e7-11eb-3e61-f96c900a0636
# ╠═87bfc240-12e3-11eb-03cc-756dc00efa6c
# ╠═c0e46442-27fb-11eb-2c94-15edbda3f84d
# ╟─bb084ace-12e2-11eb-2dfc-111e90eabfdd
# ╠═df706ebc-2a63-11eb-0b09-fd9f151cb5a8
# ╠═e2e4cfac-2a63-11eb-1b7f-9d8d5d304b43
# ╠═ecaab27e-2a16-11eb-0e99-87c91e659cf3
# ╠═e3ee80c0-12dd-11eb-110a-c336bb978c51
# ╠═9c8a7e5a-12dd-11eb-1b99-cd1d52aefa1d
# ╠═d96c7a56-12e4-11eb-123c-d57487bd37df
# ╟─6b3b6030-2066-11eb-3343-e19284638efb
