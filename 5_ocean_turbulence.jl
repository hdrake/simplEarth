### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 0bb3e89a-0d82-11eb-0d51-e7e8c38cedbb
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add("Plots")
	using Plots
end

# ╔═╡ 1a2343ee-0d82-11eb-3b6f-6fadae311344
md"""

### A simple but turbulent ocean model

Potential vorticity, the fluid-dynamic equivalent of angular momentum, is a useful quantity for oceanography. It is defined as 

$q \equiv R\zeta + y,$
where $R \equiv \pi \frac{\tau₀}{\rho H \beta^2 L^3} $ is a constant, called the Rossby number, $\zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}$ is the *relative vorticity*, and $f=y$ is the Coriolis parameter.

The vorticity equation
\begin{gather}
R^{-1} \frac{\partial q}{\partial t} + u \frac{\partial q}{\partial x} + v\frac{\partial q}{\partial y} = F - D,
\end{gather}
or 
\begin{gather}
R^{-1} \frac{\partial q}{\partial t} + J(\psi, q) = F - D,
\end{gather}
where the streamfunction $\psi$ is found by solving the Poisson equation

$\zeta = \nabla^{2}\psi,$
"""

# ╔═╡ bf7c1c06-29b8-11eb-2b87-d16032f11b13
md"""

1) Classical Munk 1966 formulation

2) Marshall 1984 formulation

"""

# ╔═╡ a586ae48-29b7-11eb-330d-cd4b9686b6e6
md"""
Boundary conditions are that

A) No-flow: $ψ = 0$ on all boundaries

B) No-slip: $\frac{\partial v}{\partial x} = 0$ (or $\frac{\partial \zeta}{\partial x} = 0$) on zonal boundaries and $\frac{\partial u}{\partial y} = 0$ ($\frac{\partial \zeta}{\partial y} = 0$) on meridional boundaries.
"""

# ╔═╡ 6c69565a-29b8-11eb-32d5-634cc9835d3b
md""" 
#### 1.2) Numerical solution method

Invert the Poisson equation $\zeta = \nabla^{2} \psi$ using the Jacobi method.

Time-step the relative vorticity using Arakawa's discretization:

Leap-frogging in time
"""

# ╔═╡ e4a27062-29bc-11eb-12c1-2577c0d3ffcd
md"""
\begin{align}
& J_{i,j}(\zeta, \psi) = - \frac{1}{12 d^2} \bigg[ & \newline
&
& \; \left( \psi _{i,j-1} + \psi _{i+1,j-1} - \psi _{i,j+1} - \psi _{i+1,j+1} \right)\left( \zeta _{i+1,j} - \zeta _{i,j} \right) \newline
&
& + \left( \psi _{i-1,j-1} + \psi _{i,j-1} - \psi _{i-1,j+1} - \psi _{i,j+1} \right)\left( \zeta _{i,j} - \zeta _{i-1,j} \right) \newline
&\bigg] &
\end{align}
"""

# ╔═╡ b304d45c-29bc-11eb-2b24-57dbaf0ae47c
md"$\left( \psi_{i,j-1} + \psi_{i+1,j-1} - \psi_{i,j+1} \right)$"

# ╔═╡ 776c2c40-0d83-11eb-1824-b3b471bd61bf
begin
	nx = 7;
	ny = 9;
	
	δt = 0.001
	δx = 0.1
	δy = 0.1
	
	Lx = δx*nx
	Ly = δy*ny
	
	x = δx*(0.5:1.:nx)'
	y = δy*(0.5:1.:ny)
end;

# ╔═╡ b86a64a6-0d8c-11eb-1cff-6bc28b7bb559
begin
	∂x(ϕ) = (ϕ[:,2:end] - ϕ[:,1:end-1])/δx
	∂y(ϕ) = (ϕ[2:end,:] - ϕ[1:end-1,:])/δy
	
	xpad(ϕ) = hcat(zeros(size(ϕ,1)), ϕ, zeros(size(ϕ,1)))
	ypad(ϕ) = vcat(zeros(size(ϕ,2))', ϕ, zeros(size(ϕ,2))')
	
	function diagnose_velocities(ψ)
		u = ∂y(ψ)
		v = ∂x(ψ)
		return u,v
	end
end

# ╔═╡ f4ec1aaa-0d8c-11eb-3024-3f1630bf9a88
function advection_stencil(q, ψ, i, j)
	-1/12*(δx*δy) * q[i+1]
end

# ╔═╡ 6449ece8-0d87-11eb-002f-b7929d86112a
function advect!(q)
	∂q∂t = zeros((nx, ny))
	
end

# ╔═╡ 2da844c0-0d85-11eb-00b6-9327f0eb2032


# ╔═╡ 9317fa80-0d83-11eb-145a-8132c1bc0bd1
Gaussian(x,y; x0=0., y0=0., σx=1., σy=1.) = exp.(-((x.-x0)/σx).^2 .-((y.-y0)/σy).^2)

# ╔═╡ e4cf21dc-0d82-11eb-0c3d-f949de09fe3a
begin
	q = Gaussian(x, y*0. .+ Ly/2., x0=Lx/2., y0=Ly/2., σx=0.2, σy=0.3)
	ψ = Gaussian(x, y, x0=Lx/2., y0=Ly/2., σx=0.3, σy=0.3)
end

# ╔═╡ 56dd1c22-0d89-11eb-3e03-b14b380dbc92
u, v = diagnose_velocities(ψ)

# ╔═╡ 44fd3fa4-0d8a-11eb-087d-3546ddd47600
size(u), size(v)

# ╔═╡ afdaf0c4-0d89-11eb-2fd6-859e33daf4c3
begin
	pu = plot(xpad(u), st=:heatmap)
	pv = plot(v, st=:heatmap)
	plot(pu, pv, layout = (1, 2), legend = false, size=(700, 280))
end

# ╔═╡ 9050a2e8-0d85-11eb-127b-eb384307258e
begin
	p1 = plot(x', y, q, st=:heatmap)
	p2 = plot(x', y, ψ, st=:heatmap)
	plot(p1, p2, layout = (1, 2), legend = false, size=(700, 280))
end

# ╔═╡ Cell order:
# ╟─1a2343ee-0d82-11eb-3b6f-6fadae311344
# ╠═bf7c1c06-29b8-11eb-2b87-d16032f11b13
# ╟─a586ae48-29b7-11eb-330d-cd4b9686b6e6
# ╠═6c69565a-29b8-11eb-32d5-634cc9835d3b
# ╠═e4a27062-29bc-11eb-12c1-2577c0d3ffcd
# ╠═b304d45c-29bc-11eb-2b24-57dbaf0ae47c
# ╠═776c2c40-0d83-11eb-1824-b3b471bd61bf
# ╠═e4cf21dc-0d82-11eb-0c3d-f949de09fe3a
# ╠═b86a64a6-0d8c-11eb-1cff-6bc28b7bb559
# ╠═56dd1c22-0d89-11eb-3e03-b14b380dbc92
# ╠═44fd3fa4-0d8a-11eb-087d-3546ddd47600
# ╠═afdaf0c4-0d89-11eb-2fd6-859e33daf4c3
# ╠═f4ec1aaa-0d8c-11eb-3024-3f1630bf9a88
# ╠═6449ece8-0d87-11eb-002f-b7929d86112a
# ╟─2da844c0-0d85-11eb-00b6-9327f0eb2032
# ╠═9050a2e8-0d85-11eb-127b-eb384307258e
# ╠═9317fa80-0d83-11eb-145a-8132c1bc0bd1
# ╠═0bb3e89a-0d82-11eb-0d51-e7e8c38cedbb
