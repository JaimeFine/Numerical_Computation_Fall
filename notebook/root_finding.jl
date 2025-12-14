### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ bd275b60-d803-11f0-839e-95ab67e9dace
md"""
# Root-Finding Algorithms

This is a notebook about root-finding algorithms in numerical computation!

## Roadmap of Root-Finding Methods

### Methods Covered

- **Bisection Method**  

- **Fixed Point Iteration**  

- **Relaxation Method**  

- **Aitken Acceleration (Steffensen’s Method)**  

- **Newton–Raphson Method**  

- **Simplified Newton Method**  

- **Secant Method**  

- **Damped Newton Method**  

---

### Comparison at a Glance

| Method                  | Needs Derivative | Needs Interval | Convergence | Robustness |
|--------------------------|------------------|----------------|-------------|------------|
| Bisection               | No               | Yes            | Linear      | Very high  |
| Fixed Point             | No               | No             | Linear      | Moderate   |
| Relaxation              | No               | No             | Linear      | Adjustable |
| Steffensen (Aitken)     | No               | No             | Quadratic   | Moderate   |
| Newton–Raphson          | Yes              | No             | Quadratic   | Sensitive  |
| Simplified Newton       | Yes (once)       | No             | Linear      | Moderate   |
| Secant                  | No               | Two guesses    | Superlinear | Moderate   |
| Damped Newton           | Yes              | No             | Quadratic   | Safer      |

This roadmap helps you navigate the trade‑offs: robustness vs speed, derivative vs derivative‑free, bracketing vs open methods.

-----

## Example

An example with function $x^{3}-x-1=0$

### Bisection Method
"""

# ╔═╡ 557ba160-8835-4c39-8af8-35c789c066e0
begin
	# Define the function
	f(x) = x^3 - x - 1
	
	function bisection(f, a, b; tol=1e-6)
	    c = (a + b) / 2 	# c is the output
	    count = 0
	
	    while (b - a) > tol
	        x = (a + b) / 2
	        t = f(x)
	        c = x
	
	        if abs(t) < tol
	            break
	        elseif t * f(a) < 0
	            b = x
	        else
	            a = x
	        end
	    end
	
	    println("Approximation: $c")
	    return c
	end
	
	bisection(f, 1.0, 2.0; tol=1e-6)
end

# ╔═╡ 8d65558a-f83d-471f-a2a8-ccbbabee95dd
md"""
### Fixed Point Iteration (Recursive Approaches)
"""

# ╔═╡ ac69754f-3d13-449b-89b5-9faff981f69d
begin
	# Define g(x) for f(x) = x^3 - x - 1
	g(x) = (x + 1)^(1/3)

	# N is max iterations L is Lipschitzz constant:
	function fixed_point(g, x0; N=100, L=0.5, tol=1e-6)
	    count = 0
	    flag = 0
	    e_refined = (1 - L) * tol
	
	    while count <= N
	        x = g(x0)
	        count += 1
	        if abs(x - x0) < e_refined
	            flag = 1
	            println("Converged at iteration $count, approximation: $x")
	            return x
	        end
	        x0 = x
	    end
	
	    println("Failed to converge within $N iterations, approximation: $x")
	    return nothing
	end
	
	fixed_point(g, 1.5; N=100, L=0.5, tol=1e-6)

end

# ╔═╡ b3d326f4-8c33-4afe-ae40-e82fa9ab4b86
md"""
### Relaxation Method*
"""

# ╔═╡ 73ca4e17-e96a-40f1-afb7-5e4fa34637b0
begin
	# g(x) was already declared previously as: g(x) = (x + 1)^(1/3)

	function relaxation(g, x0; w=1.0, tol=1e-6, N=100)
		count = 0
		while count < N
			x = (1 - w) * x0 + w * g(x0)
			count += 1
			if abs(x - x0) < tol
				println("Converged at iteration $count, approximation: $x")
				return x
			end
			x0 = x
		end
		println("Failed to converge within $N iterations, approximation: $x")
		return nothing
	end

	relaxation(g, 1.5)
end

# ╔═╡ 7635ddb9-34ed-46ba-ac2a-62addda69fd7
md"""
### Aitken Acceleration (Steffensen's Method)
"""

# ╔═╡ 4ea3a8b2-2101-4206-ae17-7e4ea79d5577
begin
	# g(x) was already declared previously as: g(x) = (x + 1)^(1/3)

	function steffensen(g, x0; tol=1e-6, N=100)
		count = 0
		while count < N
			x1 = g(x0)
			x2 = g(x1)
			# Aitken's formula here:
			denom = x2 - 2x1 + x0
			if denom == 0
				println("Denominator cannot be 0!")
			end
			
			x = x0 - (x1 - x0)^2 / denom
			count += 1
			if abs(x - x0) < tol
				println("Converged at iteration $count, approximation: $x")
				return x
			end
			x0 = x
		end
		println("Failed to converge within $N iterations, approximation: $x")
		return nothing
	end

	steffensen(g, 1.5)
end

# ╔═╡ 226b6262-cfba-4609-88e6-1f5d488ed1bd
md"""
### Newton-Raphson Method
"""

# ╔═╡ dc99680a-c025-45e3-8ea8-0c1d99324915
begin
	# f(x) was already declared previously as: f(x) = x^3 - x - 1
	# Declaration of derivative:
	df(x) = 3x^2 - 1

	function newton(f, df, x0; tol=1e-6, N=100)
		count = 0
		while count < N
			x = x0 - f(x0) / df(x0)
			if abs(x - x0) < tol
				println("Converged at iteration $count, approximation: $x")
				return x
			end
			x0 = x
			count += 1
		end
		println("Failed to converge within $N iterations, approximation: $x")
		return nothing
	end

	newton(f, df, 2)
end

# ╔═╡ b218b3d8-2070-4187-bec8-03584a58802e
md"""
### Simplified Newton Method
"""

# ╔═╡ 16fd1b49-d501-4a58-a3bd-588e2780fd8f
begin
	# f(x) and df(x) was already declared previously:
	# f(x) = x^3 - x - 1
	# df(x) = 3x^2 - 1

	function simplified_newton(f, df, x0; tol=1e-6, N=100)
		d = df(x0)
		count = 0
		while count < N
			x = x0 - f(x0) / d
			if abs(x - x0) < tol
				println("Converged at iteration $count, approximation: $x")
				return x
			end
			x0 = x
			count += 1
		end
		println("Failed to converge within $N iterations, approximation: $x")
		return nothing
	end

	simplified_newton(f, df, 2)
end

# ╔═╡ 2d5218d6-176d-440b-8f26-e17fc2a86352
md"""
### Secant Method
"""

# ╔═╡ ed749e0d-b474-41f8-8802-e59154840712
begin
	# f(x) was already declared previously as: f(x) = x^3 - x - 1

	function secant(f, x0, x1; tol=1e-6, N=100)
		count = 0
		while count < N
			if f(x1) == f(x0)
	            println("Division by zero!!!")
	            return nothing
	        end
			denom = (x1 - x0) / (f(x1) - f(x0)) 
			x = x1 - f(x1) * denom
			count += 1
			if abs(x - x1) < tol
				println("Converged at iteration $count, approximation: $x")
				return x
			end
			x0, x1 = x1, x
		end
		println("Failed to converge within $N iterations, approximation: $x")
		return nothing
	end

	secant(f, 1, 2)
end

# ╔═╡ f2f8ebc6-6975-48a6-a515-4cbb42226964
md"""
### Damped Newton
"""

# ╔═╡ f5810624-ec67-47ba-9027-39ca6017e5a6
begin
	# f(x) and df(x) was already declared previously:
	# f(x) = x^3 - x - 1
	# df(x) = 3x^2 - 1

	function newton_damped(f, df, x0; tol=1e-6, N=100)
		count = 0
		lambdas = [1 / 2^k for k in 0:10]
		
		while count < N
			for λ in lambdas
				x = x0 - λ * f(x0) / df(x0)
				if abs(x - x0) < tol
                	println("Converged at iteration $count with λ=$λ, approximation: $x")
    	            return x
        	    end
				x0 = x
			end
			count += 1
		end
		
		println("Failed to converge within $N iterations, approximation: $x")
		return nothing
	end

	newton_damped(f, df, 2)
end

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
# ╟─bd275b60-d803-11f0-839e-95ab67e9dace
# ╠═557ba160-8835-4c39-8af8-35c789c066e0
# ╟─8d65558a-f83d-471f-a2a8-ccbbabee95dd
# ╠═ac69754f-3d13-449b-89b5-9faff981f69d
# ╟─b3d326f4-8c33-4afe-ae40-e82fa9ab4b86
# ╠═73ca4e17-e96a-40f1-afb7-5e4fa34637b0
# ╟─7635ddb9-34ed-46ba-ac2a-62addda69fd7
# ╠═4ea3a8b2-2101-4206-ae17-7e4ea79d5577
# ╟─226b6262-cfba-4609-88e6-1f5d488ed1bd
# ╠═dc99680a-c025-45e3-8ea8-0c1d99324915
# ╟─b218b3d8-2070-4187-bec8-03584a58802e
# ╠═16fd1b49-d501-4a58-a3bd-588e2780fd8f
# ╟─2d5218d6-176d-440b-8f26-e17fc2a86352
# ╠═ed749e0d-b474-41f8-8802-e59154840712
# ╟─f2f8ebc6-6975-48a6-a515-4cbb42226964
# ╠═f5810624-ec67-47ba-9027-39ca6017e5a6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
