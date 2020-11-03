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

# ╔═╡ 9c8a7e5a-12dd-11eb-1b99-cd1d52aefa1d
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	using Plots
	using PlutoUI
end

# ╔═╡ 65da5b38-12dc-11eb-3505-bdaf7834afaa
begin
	Δx = 0.025
	Δy = 0.025
	Δt = 0.001
	
	κ = 0.1
	
	x = (0. +Δx/2.:Δx:1.)'
	y = (-1. +Δy/2.:Δy:1.)
end;

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

# ╔═╡ 07a48564-12e8-11eb-0eda-17b9a22cf94c
function diffuse(T)
	tend = (
		κ.*(circshift(T, (0, 1)) .- 2*T .+ circshift(T, (0, -1)))/(2Δx) .+
		κ.*(circshift(T, (1, 0)) .- 2*T .+ circshift(T, (-1, 0)))/(2Δy)
	)
	return tend
end;

# ╔═╡ a21195b2-130e-11eb-0b04-39d88d891a4c
diff_filter = [[0, 1, 0], [1, -4, 1], [0, 1, 0]]

# ╔═╡ 440fe49a-12e5-11eb-1c08-f706f5f33c84
@bind go Button("Timestep")

# ╔═╡ 1528ed7e-12e5-11eb-34cf-112d2baa7353
function temperature_heatmap(T)
	p = contourf(x', y, T, color=:bluesreds, levels=-1:0.25:1., colorbar_title="Temperature [°C]")
	plot!(clims=(-1., 1.))
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
	xψ = (0.:Δx:1.)'
	yψ = (-1:Δy:1.)
	
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
end;

# ╔═╡ 35e38798-12fe-11eb-1e89-c70a35af576b
maximum(U)

# ╔═╡ 3b4e4722-12fe-11eb-238d-17aea2c23f58
begin
	CFL_adv = maximum(V)*Δt/Δx
	CFL_diff = κ*Δt/(Δx^2)
	CFL_adv, CFL_diff
end

# ╔═╡ 1aa84632-12e3-11eb-0c5b-ebc4315d820c
function advect(T)
	tend = (
		U.*(circshift(T, (0, 1)) .- circshift(T, (0, -1)))/(2Δx) .+
		V.*(circshift(T, (1, 0)) .- circshift(T, (-1, 0)))/(2Δy)
	)
	return tend
end;

# ╔═╡ 87bfc240-12e3-11eb-03cc-756dc00efa6c
function timestep!(t, T)
	T .+= Δt*(advect(T) .+ diffuse(T))
	t .+= Δt
end;

# ╔═╡ 3b0e16a2-12e5-11eb-3130-c763c1c85182
begin
	⏩ = nothing
	go
	nT = 40
	for i = 1:nT
		timestep!(t, T)
	end
end;

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
	quiver!(p, X[4:8:end], Y[4:8:end], quiver=(U[4:8:end]./10., V[4:8:end]./10.), color=:black, alpha=0.7)
	plot!(p, xlims=(0., 1.), ylims=(-1, 1.))
	plot!(p, xlabel="longitudinal distance", ylabel="latitudinal distance")
	plot!(p, clabel="Temperature")
	p
end |> as_svg

# ╔═╡ Cell order:
# ╠═65da5b38-12dc-11eb-3505-bdaf7834afaa
# ╠═35e38798-12fe-11eb-1e89-c70a35af576b
# ╠═3b4e4722-12fe-11eb-238d-17aea2c23f58
# ╠═c4424838-12e2-11eb-25eb-058344b39c8b
# ╟─f5ae1756-12e9-11eb-1228-8f03879c154a
# ╟─f9824610-12e7-11eb-3e61-f96c900a0636
# ╠═1aa84632-12e3-11eb-0c5b-ebc4315d820c
# ╠═07a48564-12e8-11eb-0eda-17b9a22cf94c
# ╠═a21195b2-130e-11eb-0b04-39d88d891a4c
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
