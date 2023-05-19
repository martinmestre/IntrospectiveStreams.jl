### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ aa0728bb-dd20-498b-a470-831e238e1bb2
begin
	using Pkg; Pkg.activate("../")
	using IntrospectiveStreams
end

# ╔═╡ f6dc3d94-f68d-11ed-1495-c11c11bbb65c
md"# Community Atlas of Stellar Streams" 

# ╔═╡ 5c7cee94-5226-4376-853d-3f10ba846054
md"## Searching the streams at order zero
Here we show the result of making approximate cuts in proper motion (PM) and color-magnitude (CM) for the five streams using Gaia DR3. We also show the case of GD-1 with DR2 for comparison.
The PM cut is guided by the Galstream track. The filtering is done by keeping only stars that satify:
$\forall \phi_1: |\mu_1^{(t)}-\mu_1(\phi_1^{(t)}|<0.7$ and
$\forall \phi_1: |\mu_2^{(t)}-\mu_2(\phi_1^{(t)}|<0.7$
where $\mu_{1,2}$ are the PMs not corrected by the solar reflex-motion.
"

# ╔═╡ d0252793-efdb-455b-a23a-02ff8896b34d


# ╔═╡ Cell order:
# ╟─f6dc3d94-f68d-11ed-1495-c11c11bbb65c
# ╠═aa0728bb-dd20-498b-a470-831e238e1bb2
# ╠═5c7cee94-5226-4376-853d-3f10ba846054
# ╠═d0252793-efdb-455b-a23a-02ff8896b34d
