### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 1d69c848-0dc8-11eb-166a-a5fec7c97b6a
md"""
## A simple spatial model of atmosphere-ocean climate dynamics
"""

# ╔═╡ 55a4467a-0dc8-11eb-2497-89398cd58b43
md"""
Recall the zero-dimensional energy balance model of Earth's climate from Lecture **??**:

\begin{align}
\text{\color{brown}{change in heat content}} = &\;\quad \text{\color{orange}{absorbed solar radiation}} \newline
& - \text{\color{blue}{thermal cooling to space}}
\end{align}

or 

\begin{gather}
\color{brown}{C \frac{dT}{dt}}
\; \color{black}{=} \; \color{orange}{\frac{(1 - α)S}{4}}
\; \color{black}{-} \; \color{blue}{(A + BT)}
\end{gather}

where $T(t)$ was only a function of temperature and the equation was thus an **Ordinary Differential Equation (ODE)**.

In this notebook, we will expand the model to the spatial dimensions $\mathbf{x}$ (longitude, or West-East) and $\mathbf{y}$ (latitude, or South-North), and will allow ocean currents to transport heat from one grid cell ($T_{i,j}$) to a neighboring grid cell ($T_{i \pm 1,\, j}$ or $T_{i,\, j \pm 1}$).

\begin{align}
\text{\color{brown}{change in heat content}} =
&\;\quad \text{\color{orange}{absorbed solar radiation}} \newline
& - \text{\color{blue}{thermal cooling to space}} \newline
& - \text{\color{purple}{heat export by ocean currents}} 
\end{align}

"""

# ╔═╡ a9601cec-0dca-11eb-0497-7d4a7966a459
md"""
#### Adding spatial dimensions to our energy balance model

For simplicity, we will start by adding spatial dimensions to our energy balance model without any other additional complexities, such as the transport of heat from one location to another.

"""

# ╔═╡ d0483e1c-0dc9-11eb-2db3-89434516a36d
md"""

We want to solve the equation

\begin{gather}
\color{brown}{C \frac{\partial T}{\partial t}}
\; \color{black}{=} \; \color{orange}{\frac{(1 - α)S}{4}}
\; \color{black}{-} \; \color{blue}{(A + BT)}
\end{gather}

where we now allow $T = T(x,y,t)$ to be a function of $x$ (longitude), $y$ (latitude), and $t$ (time).

The discretized form of this equation,

$T_{i,\, j,\, n+1} = T_{i,\, j,\, n} + \Delta t \left[ \frac{1}{C} \left( \frac{ \left( 1-\alpha(T_{i,\, j,\, n}) \right) S}{4} - (A + BT_{i,\, j,\, n}) \right) \right]$

shows that the future temperature $T_{i,j,n+1}$ in the `(i,j)` cell only requires knowledge of the present temperature $T_{i,j,n}$ that same cell.

Thus, the discretized two-dimensional problems simple amounts to solving the zero-dimensional energy balance model in each cell of the grid.
"""

# ╔═╡ 08937658-0dcd-11eb-2ed5-535e0256641d


# ╔═╡ 5574d9cc-0dcc-11eb-2a94-cfef3a89cd40
md"""
### Adding oceanic heat transport
"""

# ╔═╡ 714d1bf4-0dc8-11eb-0fe0-3b52dc7a4d35
md"""
Consider the temperature $T(x,y,t)$ evolution equation

$\frac{\partial T}{\partial t} = - \frac{\partial}{\partial x}\left( uT - \kappa \frac{\partial T}{\partial x} \right) - \frac{\partial}{\partial y}\left( vT - \kappa \frac{\partial T}{\partial y} \right) + \frac{(1 - α)S}{4C} - \left(\frac{A}{C} + \frac{B}{C}T \right)$
"""

# ╔═╡ f89ca0ba-0dce-11eb-3327-ad9b1af679ac
md"""
##### Velocity field of a typical subtropic ocean gyre

TO DO: Code up analytical solution from Vallis or Pedloskly textbook
"""

# ╔═╡ af1ed97c-0dcd-11eb-3357-338aa6045322
# Quiver plots of velocity field (u,v)

# ╔═╡ 643ecc8a-0dcf-11eb-1b50-f30a098dfff7
md"""
##### Building intuition about heat transport
"""

# ╔═╡ 343ff15c-0dd0-11eb-0d69-7179488c11d1
# Heat maps of heat transport in x and y

# Heat maps of heat convergence in x and y

# Heat map of buoyuancy flux magnitude & visualization of divergence (?)

# ╔═╡ 0be9a322-0dcd-11eb-1258-97300c0b9d65
md"""
##### Adding a latitudinal dependence to the solar insolation

TO DO: Add a simple formula for the latitudinal dependence of insolation (in W/m$^2$/m ?)
"""

# ╔═╡ ae55d0d4-0dcd-11eb-30bc-c19b44a332fc


# ╔═╡ Cell order:
# ╟─1d69c848-0dc8-11eb-166a-a5fec7c97b6a
# ╠═55a4467a-0dc8-11eb-2497-89398cd58b43
# ╟─a9601cec-0dca-11eb-0497-7d4a7966a459
# ╟─d0483e1c-0dc9-11eb-2db3-89434516a36d
# ╠═08937658-0dcd-11eb-2ed5-535e0256641d
# ╟─5574d9cc-0dcc-11eb-2a94-cfef3a89cd40
# ╠═714d1bf4-0dc8-11eb-0fe0-3b52dc7a4d35
# ╠═f89ca0ba-0dce-11eb-3327-ad9b1af679ac
# ╠═af1ed97c-0dcd-11eb-3357-338aa6045322
# ╠═643ecc8a-0dcf-11eb-1b50-f30a098dfff7
# ╠═343ff15c-0dd0-11eb-0d69-7179488c11d1
# ╟─0be9a322-0dcd-11eb-1258-97300c0b9d65
# ╠═ae55d0d4-0dcd-11eb-30bc-c19b44a332fc
