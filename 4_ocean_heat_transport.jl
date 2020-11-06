### A Pluto.jl notebook ###
# v0.12.7

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
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	Pkg.add("ImageIO")
	Pkg.add("Images")
	Pkg.add("ImageFiltering")
	Pkg.add("OffsetArrays")
	using Plots
	using PlutoUI
	using ImageIO, Images, ImageFiltering
	using OffsetArrays
end

# ╔═╡ ed741ec6-1f75-11eb-03be-ad6284abaab8
html"""
<iframe width="700" height="394" src="https://www.youtube.com/embed/H4HUJs6LQfI" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 65da5b38-12dc-11eb-3505-bdaf7834afaa
begin
	Δx = 0.02
	Δy = 0.02
	Δt = 0.001
	
	κ = 0.1
	
	x = (0. -Δx/2.:Δx:1. +Δx/2.)'
	y = (-1. -Δy/2.:Δy:1. +Δx/2.)
	
	Nx = size(x, 2)
	Ny = size(y, 1)
end;

# ╔═╡ 9036dc6a-204e-11eb-305d-45e760e62bef
begin
	diff_kernel = OffsetArray(zeros(3,3), -1:1, -1:1)
	diff_kernel[0, 0] = -4
	diff_kernel[-1, 0] = 1; diff_kernel[1, 0] = 1;
	diff_kernel[0, -1] = 1; diff_kernel[0, 1] = 1;
	diff_kernel
end

# ╔═╡ 79a0086c-2050-11eb-1974-49d430b5eecd
begin
	function diffuse(T, j, i)
		return κ.*sum(diff_kernel[-1:1,-1:1].*T[j-1:j+1, i-1:i+1])/(2Δx^2)
	end
	diffuse(T) = [diffuse(deepcopy(T), j, i) for j=2:Ny-1, i=2:Nx-1]
end

# ╔═╡ 1cea2b90-205d-11eb-0d06-7df64faf1b53
begin
	adv_kernel = OffsetArray(zeros(3,3), -1:1, -1:1)
	adv_kernel[-1, 0] = -1; adv_kernel[1, 0] = 1;
	adv_kernel[0, -1] = -1; adv_kernel[0, 1] = 1;
	adv_kernel
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
	T = -repeat(y, 1, size(x,2));
	t = [0.]
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
@bind go Button("Timestep")

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

# ╔═╡ 989df49a-205c-11eb-25b0-4d0119d3adef
begin
	function advect(T, j, i)
		return .-(
			sum(adv_kernel[0, -1:1].*(U[j, i-1:i+1].*T[j, i-1:i+1]))/(2Δx) .+
			sum(adv_kernel[-1:1, 0].*(V[j-1:j+1, i].*T[j-1:j+1, i]))/(2Δy)
		)
	end
	advect(T) = [advect(deepcopy(T), j, i) for j=2:Ny-1, i=2:Nx-1]
end

# ╔═╡ 87bfc240-12e3-11eb-03cc-756dc00efa6c
function timestep!(t, T)
	update_ghostcells!(T)
	T[2:end-1, 2:end-1] .+= Δt*(advect(T) .+ diffuse(T))
	t .+= Δt
end;

# ╔═╡ 3b0e16a2-12e5-11eb-3130-c763c1c85182
begin
	⏩ = nothing
	go
	nT = 20
	for i = 1:nT
		timestep!(t, T)
	end
end;

# ╔═╡ 3b4e4722-12fe-11eb-238d-17aea2c23f58
begin
	CFL_adv = maximum(V)*Δt/Δx
	CFL_diff = κ*Δt/(Δx^2)
	CFL_adv, CFL_diff
end

# ╔═╡ 3cc1218e-1307-11eb-1907-e7cd68f6af35
heatmap(x', y, ψ̂)

# ╔═╡ d96c7a56-12e4-11eb-123c-d57487bd37df
as_svg(x) = PlutoUI.Show(MIME"image/svg+xml"(), repr(MIME"image/svg+xml"(), x))

# ╔═╡ bd879bbe-12de-11eb-0d1d-93bba42b6ff9
begin
	⏩
	X = repeat(xitp(x), size(yitp(y),1), 1)
	Y = repeat(yitp(y), 1, size(xitp(x),2))
	p = temperature_heatmap(T)
	Nq = 8
	quiver!(p, X[(Nq+1)÷2:Nq:end], Y[(Nq+1)÷2:Nq:end], quiver=(U[(Nq+1)÷2:Nq:end]./10., V[(Nq+1)÷2:Nq:end]./10.), color=:black, alpha=0.7)
	plot!(p, xlims=(0., 1.), ylims=(-1.0, 1.0))
	plot!(p, xlabel="longitudinal distance", ylabel="latitudinal distance")
	plot!(p, clabel="Temperature")
	p
end |> as_svg

# ╔═╡ Cell order:
# ╟─ed741ec6-1f75-11eb-03be-ad6284abaab8
# ╠═65da5b38-12dc-11eb-3505-bdaf7834afaa
# ╠═9036dc6a-204e-11eb-305d-45e760e62bef
# ╠═79a0086c-2050-11eb-1974-49d430b5eecd
# ╠═1cea2b90-205d-11eb-0d06-7df64faf1b53
# ╠═989df49a-205c-11eb-25b0-4d0119d3adef
# ╠═b68ca886-2053-11eb-2e39-35c724ed3a3c
# ╠═c4424838-12e2-11eb-25eb-058344b39c8b
# ╠═3b4e4722-12fe-11eb-238d-17aea2c23f58
# ╟─f5ae1756-12e9-11eb-1228-8f03879c154a
# ╟─f9824610-12e7-11eb-3e61-f96c900a0636
# ╠═87bfc240-12e3-11eb-03cc-756dc00efa6c
# ╟─440fe49a-12e5-11eb-1c08-f706f5f33c84
# ╠═bd879bbe-12de-11eb-0d1d-93bba42b6ff9
# ╠═3cc1218e-1307-11eb-1907-e7cd68f6af35
# ╠═3b0e16a2-12e5-11eb-3130-c763c1c85182
# ╠═1528ed7e-12e5-11eb-34cf-112d2baa7353
# ╠═bb084ace-12e2-11eb-2dfc-111e90eabfdd
# ╠═627eb1a4-12e2-11eb-30d1-c1ad292d1522
# ╠═e3ee80c0-12dd-11eb-110a-c336bb978c51
# ╠═9c8a7e5a-12dd-11eb-1b99-cd1d52aefa1d
# ╠═d96c7a56-12e4-11eb-123c-d57487bd37df
