### A Pluto.jl notebook ###
# v0.12.4

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
where $R$ is a constant, called the Rossby number, $\zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}$ is the *relative vorticity*, and $f=y$ is the Coriolis parameter.

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

Let us first consider the conservative form of this equation– no forcing and no dissipation.

"""

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
