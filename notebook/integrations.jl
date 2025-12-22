### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ e26e14f0-de15-11f0-b254-274de0af11a6
md"""
# Numerical Integrations

This is a notebook about integration algorithm in numerical computation!

## Example


### Trapezoidal Rule (Newton-Cotes, degree 1)
"""

# ╔═╡ 1367f244-51b8-4070-be00-5e7cf5269d68
begin
	function trapezoid(f, a, b, n)
		h = (b - a) / n
		total = 0.0

		for i in 0:n-1
			x0 = a + i * h
			x1 = x0 + h
			total += (h / 2) * (f(x0) + f(x1))
		end

		return total
	end

	f(x) = exp(-x^2)
	a, b = 0.0, 1.0
	n = 10

	println("Trapezoid (NC deg 1): ", trapezoid(f, a, b, n))
end

# ╔═╡ d12a4947-3dbb-4521-b826-18d8a9451af7
md"""
### Simpson's Rule (Newton-Cotes, degree 2)
"""

# ╔═╡ 51de4a86-f46e-4c9e-91cd-09e99ae03e39
begin
	function simpson(f, a, b, n)
		# Not that n must be even!
		h = (b - a) / n
		total = 0.0

		for i in 0:2:n-2
			x0 = a + i * h
			x1 = x0 + h
			x2 = x0 + 2h
			total += (h / 3) * (f(x0) + 4f(x1) + f(x2))
		end

		return total
	end

	println("Simpson (NC deg 2): ", simpson(f, a, b, n))
end

# ╔═╡ 85645321-7b5f-4a6c-946e-dbad2f150218


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.2"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╠═e26e14f0-de15-11f0-b254-274de0af11a6
# ╠═1367f244-51b8-4070-be00-5e7cf5269d68
# ╟─d12a4947-3dbb-4521-b826-18d8a9451af7
# ╠═51de4a86-f46e-4c9e-91cd-09e99ae03e39
# ╠═85645321-7b5f-4a6c-946e-dbad2f150218
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
