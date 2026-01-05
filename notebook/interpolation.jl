### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 93192e65-3b62-41d8-9bf5-94ac1a4f84ce
begin
	using Plots
	
	# Declare the points:
	x = [-1, 0, 2, 3]
	y = [-4, -1, 0, 3]

	function lagrange(x, y)
		n = length(x)

		return t -> begin
			s = 0.0
			for i in 1:n
				Li = 1.0
				for j in 1:n
					if j != i
						Li *= (t - x[j]) / (x[i] - x[j])
					end
				end
				s += y[i] * Li
			end
			s
		end
	end

	lagrange_value = lagrange(x, y)
	println("f(1.5) = ", lagrange_value(1.5))

	x_lagrange = -1.5:0.01:3.5
	y_lagrange = lagrange_value.(x_lagrange)
	p1 = plot(
	    x_lagrange, y_lagrange, label="LaGrange",
	    lw=2, color=:blue, ls=:dashdot
	)
end

# ╔═╡ 84456c50-dbd3-11f0-bcea-b987cc31182e
md"""
# Interpolations

This is a notebook about interpolation algorithm in numerical computation!

## Example

An example with values


|   |   |   |   |   |
|---|---|---|---|---|
| f(x) | -4 | -1 | 0 | 3 |
| x    | -1 | 0  | 2 | 3 |

Find f(x) where x = 1.5. Then plot the graph!

### LaGrange Interpolation
"""

# ╔═╡ 38ccfa36-aee0-4e3e-8b84-7d3a0ec42698
md"""
### Newton Interpolation
"""

# ╔═╡ 58572904-8c9a-4537-9ee5-5356ff8a74d3
begin
	# x and y are declared previously...
	function newton(x, y)
		n = length(x)
		a = float.(y)

		# Compute divided differences:
		for j in 2:n
			for i in n:-1:j
				a[i] = (a[i] - a[i-1]) / (x[i] - x[i-j+1])
			end
		end

		return t -> begin
			output = a[1]
			for i in 2:n
				term = a[i]
				for j in 1:(i-1)
					term *= (t - x[j])
				end
				output += term
			end
			output
		end

		"""
		>>> A more efficient evaluation: Horner's method:
		return t -> begin
		    output = a[n]
		    for i in (n-1):-1:1
		        output = output * (t - x[i]) + a[i]
		    end
		    output
		end
		"""
		
	end

	newton_value = newton(x, y)
	println("f(1.5) = ", newton_value(1.5))

	x_newton = -1.5:0.01:3.5
	y_newton = newton_value.(x_newton)
	p2 = plot(
		x_newton, y_newton, label="Newton",
		lw=2, color=:red, ls=:dash
	)
end

# ╔═╡ 736a9385-aeeb-454d-9ebf-8aea5f9f94a1
md"""
### Hermite Interpolation
"""

# ╔═╡ 1f25a63d-bad2-4a14-b519-dd6945659aa1
begin
	function hermite(x, y)
		n = length(x)
		z = repeat(x, inner=2)
		Q = zeros(2n, 2n)
		Q[:, 1] = repeat(y, inner=2)

		# Compute derivative first of all:
		dy = zeros(n)		
		for i in 2:n-1 
			dy[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
		end
		dy[1] = (y[2] - y[1]) / (x[2] - x[1])
		dy[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])
		
		for i in 1:n
			Q[2i-1, 2] = dy[i]
			Q[2i, 2] = dy[i]
		end

		for j in 3:2n
			for i in j:2n
				Q[i, j] = (Q[i, j-1] - Q[i-1, j-1]) / (z[i] - z[i-j+1])
			end
		end

		# Another method apart from return t-> ...
		function H(t)
			result = Q[2n, 2n]
			for i in (2n-1):-1:1
				result = result * (t - z[i]) + Q[i, i]
			end
			return result
		end

		return H
	end
	
	hermite_value = hermite(x, y)
	println("f(1.5) = ", hermite_value(1.5))

	x_hermite = -1.5:0.01:3.5
	y_hermite = hermite_value.(x_hermite)
	p3 = plot(
		x_hermite, y_hermite, label="Hermite",
		lw=2, color=:green, ls=:dashdotdot
	)
end

# ╔═╡ 3ca0b360-cab5-445b-89c9-7c83b29f6c6f
md"""
### Least Square Method
"""

# ╔═╡ db90b855-83e7-4c0c-83d1-da8f18e708a1
begin
	function least_square(x, y)
		degree = 3
		n = length(x)
		X = Array{Float64}(undef, n, degree+1)
		
		for i in 1:n
			for j in 0:degree
				X[i, j+1] = x[i]^j
			end
		end

		Xt = transpose(X)
		XtX = Xt * X
		Xty = Xt * y

		XtX_inv = inv(XtX)
		result = XtX_inv * Xty

		println("Coefficients: ", result)

		return x_val -> sum(result[j+1] * x_val^j for j in 0:degree)
	end

	ls_value = least_square(x, y)
	println("f(1.5) = ", ls_value(1.5))

	x_least_sq = -1.5:0.01:3.5
	y_least_sq = ls_value.(x_least_sq)
	p5 = plot(
		x_least_sq, y_least_sq, label="Least Square",
		lw=2, color=:violet, ls=:solid
	)
end

# ╔═╡ 42727f41-5c33-4bf2-9351-0f4a1f193a9c

# ╔═╡ Cell order:
# ╠═84456c50-dbd3-11f0-bcea-b987cc31182e
# ╠═93192e65-3b62-41d8-9bf5-94ac1a4f84ce
# ╟─38ccfa36-aee0-4e3e-8b84-7d3a0ec42698
# ╠═58572904-8c9a-4537-9ee5-5356ff8a74d3
# ╟─736a9385-aeeb-454d-9ebf-8aea5f9f94a1
# ╠═1f25a63d-bad2-4a14-b519-dd6945659aa1
# ╠═3ca0b360-cab5-445b-89c9-7c83b29f6c6f
# ╠═db90b855-83e7-4c0c-83d1-da8f18e708a1
# ╠═42727f41-5c33-4bf2-9351-0f4a1f193a9c
