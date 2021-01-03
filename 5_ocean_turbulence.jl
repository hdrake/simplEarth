### A Pluto.jl notebook ###
# v0.12.15

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

# ╔═╡ 0bb3e89a-0d82-11eb-0d51-e7e8c38cedbb
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		"Plots",
		"PlutoUI",
		"OffsetArrays"
	])
	using Plots
	using PlutoUI
	using OffsetArrays
end

# ╔═╡ ca5dd220-33d3-11eb-1bd1-3da9995278b6
md"""

# MIT 18.S191 – Lecture 24

## Building a simple non-linear ocean model

**Part V of the Climate Modelling unit for MIT 18.S191**

*Guest lecturer:* Henri Drake (PhD student in MIT Course 12)
"""

# ╔═╡ bfe3b0c6-33d3-11eb-2f75-499f5eafb161
html"""
<iframe width="680" height="430" src="https://www.youtube-nocookie.com/embed/CRzoCkSJFX0?start=109" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 1a2343ee-0d82-11eb-3b6f-6fadae311344
md"""

### 0. Background

##### 0.1 Review

So far we have only considered the thermodynamic (or heat) equation that tells us the temperature of a fluid based on its transport and thermal forcings, for example:

$\frac{\partial T}{\partial t} = -u\frac{\partial T}{\partial x} -v\frac{\partial T}{\partial y} + \kappa \left( \frac{\partial^{2} T}{\partial x^{2}} + \frac{\partial^{2} T}{\partial y^{2}} \right) + \mathcal{F},$

where $\mathbf{u} = u \hat{\mathbf{x}} + v \hat{\mathbf{y}}$ is a given velocity field.

But where does this velocity $\mathbf{u}$, ocean currents and atmospheric winds, actually come from?

##### 0.2 Non-linear effects (Burger's equation)

Consider the following equation

$\frac{\partial u}{\partial t} = - u \frac{\partial u}{\partial x} + \nu \frac{\partial^{2} u}{\partial x^{2}}.$

Note that this is identical to the **one-dimensional advection-diffusion** from before, but now the variable we are solving for *is the velocity* $u$ *itself*!

Using the product rule, we can rewrite this as

$\frac{\partial u}{\partial t} = -\frac{1}{2}\frac{\partial (u^{2})}{\partial x} + \nu \frac{\partial^{2} u}{\partial x^{2}}.$

Note that the advective term depends on a quadratic (or **non-linear**) power of the velocity $u$! This non-linearity brings a lot of interesting richness to the problem, but it also comes with a lot of complications in terms of numerical implementations.

"""

# ╔═╡ 88f8a238-33ce-11eb-3e95-ffdf1ee6fee1
begin
	iplot_slider = @bind iplot Slider(1:5:400, default=1)
	md"""
	*Drag this slider to advance forward in time* $(iplot_slider)
	"""
end

# ╔═╡ bd7b9070-33cd-11eb-2946-5923c86830b2
md"""
##### 0.3 The Navier-Stokes equations

In actuality, we *do know* the equations that govern fluid flow $\mathbf{u} = (u,v,w)$. They are the **Navier-Stokes** equations and are given by

$\rho \left(\frac{\partial \mathbf{u}}{\partial t} +\mathbf{u} \cdot \nabla  \mathbf{u} \right) =  - \nabla p + \rho g\,\hat{\mathbf{z}} + \nu \nabla^{2} \mathbf{u},$

where $p$ is the pressure, $\rho$ is the fluid density, $g$ is the gravitational acceleration, and $\nu$ is a kinematic viscosity (just like the diffusivity, but for momentum instead of heat). This equation is the equivalent of

$m \mathbf{a} = \mathbf{F}$

for a fluid, where

$\text{density} \times \textbf{acceleration} = \textbf{pressure forces} + \textbf{gravity} + \textbf{friction}$

Note that, just like in Burgers' equation above, we have a pesky non-linear term, $\mathbf{u} \cdot \nabla \mathbf{u}$, which we have to deal with. This non-linear term is what makes fluid dynamics (and climate dynamics) interesting.

For incompressible fluids like the ocean, the Navier-Stokes equations are by the **continuity** (conservation of mass) equation

$\nabla \cdot \mathbf{u} = 0$

##### 0.4 The two-dimensional vorticity equation

While one could solve the Navier-Stokes equation by **direct numerical simulation (DNS)**, this is only possible for very small problems (domains less than $100$ m $\times$ $100$ m $\times$ $100$ m large, at most) and is not helpfulf or modelling global-scale ocean currents.

The Navier-Stokes equations can be transformed into two coupled PDEs that describe depth-integrated ocean currents as follows (see course 12.800 and follow-on courses for the detailed derivation):

$\frac{\partial \zeta}{\partial t} = -\frac{\partial \psi}{\partial x}
+
\delta_{M}^{3} \left( \frac{\partial^{2} \zeta}{\partial x^{2}} + \frac{\partial^{2} \zeta}{\partial y^{2}} \right)
-
\delta_{I}^{2} \left( \frac{\partial \psi}{\partial x}\frac{\partial \zeta}{\partial y} - \frac{\partial \psi}{\partial y}\frac{\partial \zeta}{\partial x} \right)$

and 

$\zeta = \frac{\partial^{2} \psi}{\partial x^{2}} + \frac{\partial^{2} \psi}{\partial y^{2}}$

We have introduced two new variables: the **relative vorticity** $\zeta$ (or fluid angular momentum) and the **streamfunction** $\psi$, which uniquely described an incompressible flow field $(u, v) = (\psi_{y}, -\psi_{x})$.

"""

# ╔═╡ b21c6074-333c-11eb-2d34-4d693945eedb
md"""
##### 0.5 A two-step time-stepping algorithm

Unlike all of our previous examples, our problem is now defined by *two* coupled PDEs. Since the second equation does not include any time-derivative terms $\frac{\partial}{\partial t}$, it holds at any time and we call it a **diagnostic** equation.

Our solution algorithm is thus:

1. Start from an inition condition $\zeta^{0} \equiv \zeta(t=0)$ at $n=0$.
2. Invert the Poisson equation $\zeta^{n} = \nabla^{2} \psi^{n}$ to determine $\psi^{n}$ (use $\psi^{n-1}$ as a guess, if available!).
3. Timestep the vorticity equation forward to determine $\zeta^{n+1}$ from $\zeta^{n}$ and $\psi^{n}$.
4. Go back to step 2. and repeat to find $\zeta^{n+2}$ from $\zeta^{n+1}$ and $\psi^{n+1}$.


"""

# ╔═╡ 2966d15a-33d5-11eb-27eb-654fca8f64aa
md"""
### 1. Numerical implementation

##### 1.1 Setting up the discrete grid and parameters
"""

# ╔═╡ 1982dc74-333e-11eb-2f05-8318f56879b4
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

# ╔═╡ afe23b7e-337f-11eb-29f5-7505363c27ca
begin
	N = 30
	G = Grid(N, 1.)
	
	Δx = G.Δx
	Δt = 0.01
	
	δM = 0.05
	Re = 4.0
	δI = δM*Re^(1/3)
end;

# ╔═╡ 9b45d716-3343-11eb-1d07-f70447ef3ccf
md"""
##### 1.2 Inverting the Poisson equation for $\psi$

Below, we use Jacobi iteration (introduced by Prof. Edelman and John Urschel earlier in the class) to invert

$\zeta = \frac{\partial^{2} \psi}{\partial x^{2}} + \frac{\partial^{2} \psi}{\partial y^{2}}$

for $\psi$.

"""

# ╔═╡ 131767f2-3383-11eb-09a5-2db0723b3723
md"**Example: the Rankine vortex**

Throughout, we will use the Rankine vortex flow field to illustrate the relationship between vorticity $\zeta$, the streamfunction $\psi$, and the velocity field $(u,v)$."

# ╔═╡ 6cbe6022-33d7-11eb-167c-73d037e10807
md"**Let's test our Jacobi solver by inverting for inverting $\zeta$ for $\psi$ and then re-diagnosing $\zeta$ to see if we retrieve our input.**"

# ╔═╡ 179e9414-33d6-11eb-2f2d-f7de7484078e
begin
	Jacobi_slider = @bind jacbobi_plot Slider(1:20:1000, default=1)
	md"""
	*Drag this slider to advance forward in Jacobi iterations* $(Jacobi_slider)
	"""
end

# ╔═╡ 7ef4755a-33d8-11eb-3bbb-c5b0bf76462b
md"""
##### 1.3 Boundary conditions
"""

# ╔═╡ 38856ca6-3355-11eb-15d4-9ffc06486c48
md"""
We will use boundary conditions:
- Slip ($\zeta = 0$) and
- Impermeability ($\psi = 0$)
on all four boundaries, for simplicity.
"""

# ╔═╡ 90f8e9e8-33d8-11eb-3b2e-0913eb06b8f7
md"### 1.4 Tendency terms

We now implement each of the tendency terms in the vorticity equation.

##### 1.4.1 Advection of planetary vorticity"

# ╔═╡ f5a53b40-3355-11eb-0564-417584508f9c
md"""
$-\frac{\partial \psi}{\partial x} \approx \frac{\psi_{j, i+1} - \psi_{j,i-1}}{2\Delta x}$
"""

# ╔═╡ 9efa0934-33d8-11eb-03ae-7323d5ede5f3
md"##### 1.4.2 Vorticity diffusion"

# ╔═╡ 3b9db190-3356-11eb-3d80-a599053db795
md"""
$\delta_{M}^{3} \left( \frac{\partial^{2} \zeta}{\partial x^{2}} + \frac{\partial \zeta^{2}}{\partial y^{2}} \right) \approx \delta_{M}^{3} \left( \frac{\zeta_{j,i+1} + \zeta_{j,i-1} + \zeta_{j+1,i} + \zeta_{j-1,i} - 4 \zeta_{i,j}}{(\Delta x)^{2}}\right)$
"""

# ╔═╡ d525c7e8-3357-11eb-3b34-2543bc3b7de0
md"""
##### 1.4.3 Non-linear vorticity advection

We would like to be able to discretize the non-linear advection by applying centered finite differences, as we did for the other two terms

$- \delta_{I}^{2} \left( \frac{\partial \psi}{\partial x}\frac{\partial \zeta}{\partial y} - \frac{\partial \psi}{\partial y}\frac{\partial \zeta}{\partial x} \right) \approx$

$- \delta_{I}^{2} \left[ \left(\frac{\psi_{j,i+1} - \psi_{j,i-1}}{2\Delta x}\right) \left(\frac{\zeta_{j+1,i} - \zeta_{j-1,i}}{2\Delta x}\right) + \left(\frac{\psi_{j+1,i} - \psi_{j-1,i}}{2\Delta x}\right) \left(\frac{\zeta_{j,i+1} - \zeta_{j,i-1}}{2\Delta x}\right) \right]$

Unfortunately, this choice of a discretization does not conserve energy or angular momentum variance (*enstrophy*), two key quantities for fluid flows. It turns out that one specific linear combination of several possible discretizations does actually conserve energy and enstrophy,


$- \delta_{I}^{2} \left( \frac{\partial \psi}{\partial x}\frac{\partial \zeta}{\partial y} - \frac{\partial \psi}{\partial y}\frac{\partial \zeta}{\partial x} \right) \approx \text{a huge mess}$

but it is too long to write it out here (see the nightmarish julia implementaiton below).


"""

# ╔═╡ 4449eabe-3356-11eb-1150-b5eaf35bd047
begin
	function nonlinear_vorticity_advection(ψ, ζ, j, i)
		return -δI^2/(12*Δx^2)*(
			 (ψ[j-1,i  ]+ψ[j-1,i+1]-ψ[j+1,i  ]-ψ[j+1,i+1])*(ζ[j  ,i+1]-ζ[j  ,i  ])
			+(ψ[j-1,i-1]+ψ[j-1,i  ]-ψ[j+1,i-1]-ψ[j+1,i  ])*(ζ[j  ,i  ]-ζ[j  ,i-1])
			+(ψ[j  ,i+1]+ψ[j+1,i+1]-ψ[j  ,i-1]-ψ[j+1,i-1])*(ζ[j+1,i  ]-ζ[j  ,i  ])
			+(ψ[j-1,i+1]+ψ[j  ,i+1]-ψ[j-1,i-1]-ψ[j  ,i-1])*(ζ[j  ,i  ]-ζ[j-1,i  ])

			+(ψ[j  ,i+1]-ψ[j+1,i  ])*(ζ[j+1,i+1]-ζ[j  ,i  ])
			+(ψ[j-1,i  ]-ψ[j  ,i-1])*(ζ[j  ,i  ]-ζ[j-1,i-1])
			+(ψ[j+1,i  ]-ψ[j  ,i-1])*(ζ[j+1,i-1]-ζ[j  ,i  ])
			+(ψ[j  ,i+1]-ψ[j-1,i  ])*(ζ[j  ,i  ]-ζ[j-1,i+1])
		)
	end
	
	function nonlinear_vorticity_advection(ψ, ζ)
		return [
			nonlinear_vorticity_advection(ψ, ζ, j, i)
			for j=2:size(ψ,1)-1, i=2:size(ψ,2)-1
		]
	end
end;



# ╔═╡ b9ef3226-33d8-11eb-3ec1-f7821a9b826f
md"##### 1.4.4 Wind forcing"

# ╔═╡ 8f1f196e-3374-11eb-3f69-97fcddfd26df
wind_forcing() = -sin.(π*G.y[2:end-1, :]);

# ╔═╡ f2b20156-33d8-11eb-3bf9-ab20146178f7
plot(wind_forcing()[:], G.y[2:end-1], ylabel="y", size=(240, 300), xlabel="wind forcing", label=nothing)

# ╔═╡ 1d9c165a-3379-11eb-0aff-4bca060d1a92
@bind go Clock()

# ╔═╡ 6d7ebbb8-33d9-11eb-2be0-41feb3bd6f32
md"### Helper functions"

# ╔═╡ ea1a04c6-3379-11eb-2e44-459031d1fd2a
begin
	X = repeat(G.x, G.Ny, 1)
	Y = repeat(G.y, 1, G.Nx)
end;

# ╔═╡ 720a201e-3368-11eb-048c-51f3d477bfbe
begin
	diff_kernel = OffsetArray(
		[0  1  0
		 1 -4  0
		 0  1  0],
		-1:1, -1:1
	)
end;

# ╔═╡ b0467d58-3368-11eb-2a24-13226257b24e
begin 
	function vorticity_diffusion(ζ, j, i)
		return δM^3 * sum(diff_kernel[-1:1, -1:1].*ζ[j-1:j+1, i-1:i+1])/(Δx^2)
	end
	
	function vorticity_diffusion(ζ)
		return [vorticity_diffusion(ζ, j, i) for j=2:size(ζ,1)-1, i=2:size(ζ,2)-1]
	end
end;


# ╔═╡ e01e9b90-3369-11eb-1d52-63810ff3daba
begin
	xgrad_kernel = OffsetArray(reshape([-1., 0., 1.], 1, 3),  0:0, -1:1)
end;

# ╔═╡ f3809f44-3355-11eb-37d9-398d2abd0797
begin
	function planetary_vorticity_advection(ψ, j, i)
		return -sum(xgrad_kernel[0:0, -1:1].*ψ[j:j, i-1:i+1])/(2*Δx)
	end;
	
	function planetary_vorticity_advection(ψ)
		return [planetary_vorticity_advection(ψ, j, i)
			    for j=2:size(ψ,1)-1, i=2:size(ψ,2)-1]
	end
end;

# ╔═╡ 27d1f9b4-33d3-11eb-3f4e-55fe30a4776c
begin
	import Base.zeros
	zeros(G::Grid) = zeros(G.Ny, G.Nx)
end

# ╔═╡ ef92ca9a-333e-11eb-2250-83fc59ba2004
begin
	function jacobi_step(ψ::AbstractMatrix, ζ::AbstractMatrix, G::Grid)
	
		ψnew = zeros(size(ψ))
		ψnew[2:end-1, 2:end-1] = [
			((ψ[j+1, i] + ψ[j-1, i] + ψ[j, i-1] + ψ[j, i+1]) / 4
				- ζ[j, i]*(G.Δx^2)/4) # We only added this line
			for j=2:size(ψ,1)-1, i=2:size(ψ,2)-1
		]
		
		return ψnew
	end
	
	function poisson_solve!(ψ, ζ, G; ϵ=1e-6, num_steps=500, return_results=false)
	
		if return_results
			results = []
		end
			
	    for i in 1:num_steps
	
	        ψtmp = jacobi_step(ψ, ζ, G)
			
			if return_results
				push!(results, deepcopy(ψtmp))
			end
			
			if maximum(abs.(ψtmp - ψ)) < ϵ
				ψ[:,:] = ψtmp
				break
			end
	
	        ψ[:,:] = ψtmp
	    end
		
		if return_results
			return results
		end
	end
end;

# ╔═╡ 5854b0b4-336c-11eb-1f31-a1d06414bb7e
function timestep!(ζ, ψ)
	poisson_solve!(ψ, ζ, G)
	
	ζ[2:end-1, 2:end-1] .+= Δt.*(
		planetary_vorticity_advection(ψ) .+
		vorticity_diffusion(ζ) .+
		nonlinear_vorticity_advection(ψ, ζ) .+
		wind_forcing()
	)
end

# ╔═╡ d826a05e-33cd-11eb-1564-f7b429d85fea
let
	nx = 100
	Lx = 1.
	Δx = Lx/nx
	Δt = 0.001
	U = 1.	
	κ = 0.05


	x = Δx/2.:Δx:Lx
	

	u = U*sin.(2π*x);
	t = [0.]

	function advect(u)
		return 0.5*(circshift(u.^2, (1)) .- circshift(u.^2, (-1)))/(2Δx)
	end
	
	function diffuse(u)
		return κ*(circshift(u, (1)) .- 2*u .+ circshift(u, (-1)))/(Δx^2)
	end

	function timestep!(t, u)
		u .+= Δt*(advect(u))
		t .+= Δt
	end

	nT = 400
	results = []
	times = []
	for i = 1:nT
		timestep!(t, u)
		push!(results, deepcopy(u))
		push!(times, deepcopy(t[1]))
	end
	
	p = plot(x, results[iplot], ylims=(-1.2, 1.2), label=nothing)
	plot!(p, xlabel="x", ylabel="u", title="Burger's equation")
	annotate!(0.1, 1.11, text(string("t = ", round(times[iplot], digits=4))))
	p
end

# ╔═╡ b87ed172-33d7-11eb-0c9f-c9e753b8c158
begin
	ζ0 = zeros(G)
	x0, y0 = 0.5, 0.0
	idx = ((G.x .- x0).^2 .+ (G.y .- y0).^2).^0.5 .< 0.25
	ζ0[idx] .= 1
end;

# ╔═╡ ac424bb4-3341-11eb-3895-5fce079086c5
function diagnoseζ(ψ, G)
	ζverify = zeros(G)
	for j=2:G.Ny-1, i=2:G.Nx-1
		ζverify[j,i] = (
			ψ[j+1, i] +
			ψ[j-1, i] +
			ψ[j, i+1] +
			ψ[j, i-1] -
			4*ψ[j, i]
		)/G.Δx^2
	end
	
	return ζverify
end

# ╔═╡ 0f2c1aba-336b-11eb-0d38-b1f296715cfc
begin
	ψ0 = zeros(size(ζ0)) # initial guess is just zeros
	@elapsed results = poisson_solve!(ψ0, ζ0, G, ϵ=1e-8, num_steps=1000, return_results=true)
	ζ0verify = diagnoseζ(ψ0, G)
end;

# ╔═╡ 80e5682e-3379-11eb-0e52-ed2888849b79
begin
	ψ = deepcopy(zeros(G))
	ζ = deepcopy(zeros(G)) .+ 1.e-1*(rand(Float64, size(ψ)).-0.5)
	t = [0]
end;

# ╔═╡ b86a64a6-0d8c-11eb-1cff-6bc28b7bb559
begin
	∂x(ϕ) = (ϕ[:,2:end] - ϕ[:,1:end-1])/Δx
	∂y(ϕ) = (ϕ[2:end,:] - ϕ[1:end-1,:])/Δx
	
	xpad(ϕ) = hcat(zeros(size(ϕ,1)), ϕ, zeros(size(ϕ,1)))
	ypad(ϕ) = vcat(zeros(size(ϕ,2))', ϕ, zeros(size(ϕ,2))')
	
	function diagnose_velocities(ψ)
		u = ∂y(ψ)
		v = ∂x(ψ)
		return u,v
	end
end

# ╔═╡ 27cba404-3344-11eb-1c56-0f1b4f765b67
md"### Pluto magic"

# ╔═╡ 6db0c1a0-3341-11eb-350a-f70869b72128
as_svg(x) = PlutoUI.Show(MIME"image/svg+xml"(), repr(MIME"image/svg+xml"(), x))

# ╔═╡ 599d565e-333e-11eb-27f1-1177a94e9bf5
begin
	p1 = heatmap(G.x[:], G.y[:], ζ0,
		title="Initial vorticity ζ₀", clims=(0, 1.))
	p2 = heatmap(G.x[:], G.y[:], diagnoseζ(results[jacbobi_plot], G),
		title="vorticity from Poisson solve", clims=(0, 1.))
	p3 = heatmap(G.x[:], G.y[:], results[jacbobi_plot],
		title="streamfunction ψ", clims=(-0.04, 0.04), color=:bluesreds)
	p = plot(
		p1,
		p2,
		p3,
		layout=(1, 3),
		size=(900,330),
	)
	annotate!(p2, 0.25, 0.95, text(string("i = ", jacbobi_plot), :white))
	p
end |> as_svg

# ╔═╡ 56aa75d0-33d9-11eb-0b58-dd1d21191129
heatmap(planetary_vorticity_advection(ψ0), size=(260, 380)) |> as_svg

# ╔═╡ 4f943864-33d9-11eb-0955-b38962ded6b4
heatmap(vorticity_diffusion(ζ0), size=(260, 380)) |> as_svg

# ╔═╡ 29fc0314-33d9-11eb-0760-6953b92cfe98
heatmap(nonlinear_vorticity_advection(ψ0, ζ0), size=(260, 380)) |> as_svg

# ╔═╡ f9ecfe38-336c-11eb-305a-65dfa9892a4c
let	
	go
	for i=1:100
		timestep!(ζ, ψ)
		t[1] += 1
	end
	
	u, v = diagnose_velocities(ψ);
	
	plot(
		heatmap(G.x[:], G.y[:], ζ, title="Vorticity"),
		heatmap(G.x[:], G.y[:], ψ, title="Streamfunction"),
		size=(600,430),
	)
	annotate!(0.02, 0.97, text(
			string("t = ", round(t[1]*Δt)),
			color=:white, :left, 8
		))
end |> as_svg

# ╔═╡ Cell order:
# ╟─ca5dd220-33d3-11eb-1bd1-3da9995278b6
# ╟─bfe3b0c6-33d3-11eb-2f75-499f5eafb161
# ╟─1a2343ee-0d82-11eb-3b6f-6fadae311344
# ╟─88f8a238-33ce-11eb-3e95-ffdf1ee6fee1
# ╟─d826a05e-33cd-11eb-1564-f7b429d85fea
# ╟─bd7b9070-33cd-11eb-2946-5923c86830b2
# ╟─b21c6074-333c-11eb-2d34-4d693945eedb
# ╟─2966d15a-33d5-11eb-27eb-654fca8f64aa
# ╠═1982dc74-333e-11eb-2f05-8318f56879b4
# ╠═afe23b7e-337f-11eb-29f5-7505363c27ca
# ╟─9b45d716-3343-11eb-1d07-f70447ef3ccf
# ╠═ef92ca9a-333e-11eb-2250-83fc59ba2004
# ╟─131767f2-3383-11eb-09a5-2db0723b3723
# ╠═b87ed172-33d7-11eb-0c9f-c9e753b8c158
# ╟─6cbe6022-33d7-11eb-167c-73d037e10807
# ╠═0f2c1aba-336b-11eb-0d38-b1f296715cfc
# ╟─179e9414-33d6-11eb-2f2d-f7de7484078e
# ╟─599d565e-333e-11eb-27f1-1177a94e9bf5
# ╟─ac424bb4-3341-11eb-3895-5fce079086c5
# ╟─7ef4755a-33d8-11eb-3bbb-c5b0bf76462b
# ╟─38856ca6-3355-11eb-15d4-9ffc06486c48
# ╟─90f8e9e8-33d8-11eb-3b2e-0913eb06b8f7
# ╟─f5a53b40-3355-11eb-0564-417584508f9c
# ╠═f3809f44-3355-11eb-37d9-398d2abd0797
# ╟─56aa75d0-33d9-11eb-0b58-dd1d21191129
# ╟─9efa0934-33d8-11eb-03ae-7323d5ede5f3
# ╟─3b9db190-3356-11eb-3d80-a599053db795
# ╠═b0467d58-3368-11eb-2a24-13226257b24e
# ╟─4f943864-33d9-11eb-0955-b38962ded6b4
# ╟─d525c7e8-3357-11eb-3b34-2543bc3b7de0
# ╠═4449eabe-3356-11eb-1150-b5eaf35bd047
# ╟─29fc0314-33d9-11eb-0760-6953b92cfe98
# ╟─b9ef3226-33d8-11eb-3ec1-f7821a9b826f
# ╠═8f1f196e-3374-11eb-3f69-97fcddfd26df
# ╟─f2b20156-33d8-11eb-3bf9-ab20146178f7
# ╠═5854b0b4-336c-11eb-1f31-a1d06414bb7e
# ╠═1d9c165a-3379-11eb-0aff-4bca060d1a92
# ╠═80e5682e-3379-11eb-0e52-ed2888849b79
# ╟─f9ecfe38-336c-11eb-305a-65dfa9892a4c
# ╟─6d7ebbb8-33d9-11eb-2be0-41feb3bd6f32
# ╠═ea1a04c6-3379-11eb-2e44-459031d1fd2a
# ╠═b86a64a6-0d8c-11eb-1cff-6bc28b7bb559
# ╠═720a201e-3368-11eb-048c-51f3d477bfbe
# ╠═e01e9b90-3369-11eb-1d52-63810ff3daba
# ╠═27d1f9b4-33d3-11eb-3f4e-55fe30a4776c
# ╟─27cba404-3344-11eb-1c56-0f1b4f765b67
# ╠═0bb3e89a-0d82-11eb-0d51-e7e8c38cedbb
# ╠═6db0c1a0-3341-11eb-350a-f70869b72128
