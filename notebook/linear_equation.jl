### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 61958f32-a8b6-40d6-be2e-858574ab7708
begin
	using LinearAlgebra

	# Declare matrix A and vector b, Ax = b
	A = [10.0 -1.0 -2.0;
		 -1.0 10.0 -2.0;
		 -1.0 -1.0 5.0]
	
	b = [7.2, 8.3, 4.2]
	
	function jacobi(
		A::Matrix{Float64}, b::Vector{Float64};
		x=zeros(size(b)), tol=1e-8, max_iter=100
	)
		n = length(b)
		D = Diagonal(diag(A))
		
		# The following can be written as:
		# L = tril(A, -1)
		# U = triu(A, 1)
		# But just for knowing how we get it...
		L = LowerTriangular(A) - D
		U = UpperTriangular(A) - D
	
		D_inv = inv(D)
	
		for k in 1:max_iter
			x_new = -D_inv * (L + U) * x + D_inv * b
			if norm(x_new - x, Inf) < tol
				println("Approximation: ", x_new)
				return x_new
			end
			x = x_new
		end
		
		return x
	end

	jacobi(A, b)
end

# ╔═╡ 6ca1508e-d9b3-11f0-85e3-a9cca3c67ff7
md"""
# Numerical methods for Linear Equations

This is a notebook about numerical linear algebra algorithms in numerical computation!

## Example

An example with function

$\begin{cases}
10x_1 - x_2 - 2x_3 = 7.2 \\
-x_1 + 10x_2 - 2x_3 = 8.3 \\
-x_1 - x_2 + 5x_3 = 4.2
\end{cases}$

### Jacobi Method
"""

# ╔═╡ f880da5a-0dd1-4592-bf18-3d830441abd4
md"""
### Gauss-Seidel Method
"""

# ╔═╡ 41c19325-2ac9-42f7-862e-37abd505bfc8
begin
	function gauss_seidel(
		A::Matrix{Float64}, b::Vector{Float64};
		x=zeros(size(b)), tol=1e-8, max_iter=100
	)
		n = length(b)
		D = Diagonal(diag(A))
		L = LowerTriangular(A) - D
		U = UpperTriangular(A) - D

		DL_inv = inv(D + L)

		for k in 1:max_iter
			x_new = DL_inv * (b - U * x)
			if norm(x_new - x, Inf) < tol
				println("Approximation: ", x_new)
				return x_new
			end
			x = x_new
		end

		return x
	end

	# A and b was already declared previously
	gauss_seidel(A, b)
end

# ╔═╡ ce647298-2958-48fa-9361-3f48955a098f
md"""
### Successive Over Relaxation (SOR)
"""

# ╔═╡ 317227db-1cf5-4329-8441-9ad00a9b1013
begin
	function sor(
		A::Matrix{Float64}, b::Vector{Float64}, ω;
		x=zeros(size(b)), tol=1e-8, max_iter=100
	)
		n = length(b)
		I = ones(size(b))
		D = Diagonal(diag(A))
		L = LowerTriangular(A) - D
		U = UpperTriangular(A) - D

		M = D + ω * L
		N = (1 - ω) * D - ω * U
		B = inv(M) * N
		g = ω * inv(M) * b

		for k in 1:max_iter
			x_new = B * x + g 
			if norm(x_new - x, Inf) < tol
				println("Approximation: ", x_new)
				return x_new
			end
			x = x_new
		end

		return x
	end

	# Define the relaxation parameter:
	ω = 1.1
	
	# A and b was already declared previously...
	sor(A, b, ω)
end

# ╔═╡ 85485111-3749-4459-9b8b-d6b3186435e5
md"""
### Gaussian Elimination Method
"""

# ╔═╡ 16e9814d-42df-4fa8-8b44-e81565e3d335
begin
	function gaussian_elimination(
		A::Matrix{Float64}, b::Vector{Float64}
	)
		n = size(A, 1)
		Ab = hcat(A, b)		# Now we get the augmented matrix

		for i in 1:n-1
			for j in i+1:n
				factor = Ab[j, i] / Ab[i, i]
				Ab[j, i:end] .-= factor .* Ab[i, i:end]
			end
		end

		x = zeros(Float64, n)
		for i in n:-1:1
			# This part: dot(Ab[i, i+1:n], x[i+1:n]),
			# computes the sum directly!
			x[i] = (Ab[i, end] - dot(Ab[i, i+1:n], x[i+1:n])) / Ab[i, i]
		end

		println("Approximation: ", x)
		
		return x
	end

	gaussian_elimination(A, b)
end

# ╔═╡ 7aea5bfd-3757-48fe-adc7-4860e3bc68a8
md"""
### Partial Pivoting Gauss Elimination
"""

# ╔═╡ beddc0f4-18c0-4b22-aad8-d173023dc259
begin
	function partial_pivoting_ge(
		A::Matrix{Float64}, b::Vector{Float64}
	)
		n = size(A, 1)
		Ab = hcat(A, b)
		
		for i in 1:n-1
			pivot_row = argmax(abs.(Ab[i:end, i])) + i - 1
			if pivot_row != i
				Ab[i, :], Ab[pivot_row, :] = Ab[pivot_row, :], Ab[i, :]
			end

			for j in i+1:n
				factor = Ab[j, i] / Ab[i, i]
				Ab[j, i:end] .-= factor .* Ab[i, i:end]
			end
		end

		x = zeros(Float64, n)
		for i in n:-1:1
			x[i] = (Ab[i, end] - dot(Ab[i, i+1:n], x[i+1:n])) / Ab[i, i]
		end

		println("Approximation: ", x)
		
		return x
	end

	partial_pivoting_ge(A, b)
end

# ╔═╡ cb4087f2-5021-4f48-b244-1475303bfd2d
md"""
### Doolittle's LU Decomposition
"""

# ╔═╡ 2e7b26e0-c6be-46f0-bcd2-2cd092620217
begin
	function LU_decomposition(
		A::Matrix{Float64}, b::Vector{Float64}
	)
		n = size(A, 1)
		U = zeros(Float64, n, n)
		L = zeros(Float64, n, n)

		for r in 1:n
			L[r, r] = 1.0
		end
			
		U[1, :] = A[1, :]
		L[2:end, 1] = A[2:end, 1] ./ U[1, 1]

		for r in 2:n
			for j in r:n
				U[r, j] = A[r, j] - dot(L[r, 1:r-1], U[1:r-1, j])
			end
			for i in r+1:n
				L[i, r] = (A[i, r] - dot(L[i, 1:r-1], U[1:r-1, r])) / U[r, r]
			end
		end

		y = zeros(Float64, n)
		for i in 1:n
		    y[i] = (b[i] - dot(L[i, 1:i-1], y[1:i-1])) / L[i,i]
		end

		x = zeros(Float64, n)
		for i in n:-1:1
			x[i] = (y[i] - dot(U[i, i+1:n], x[i+1:n])) / U[i, i]
		end

		println(x)
		
		return nothing
	end

	LU_decomposition(A, b)
end

# ╔═╡ 09cf821b-991e-4e5c-823b-726a25414d9c
md"""
### Cholesky's LU Decomposition
"""

# ╔═╡ f97c4dab-dba6-4038-bd18-64c959a2a8e2
begin
	function cholesky(
		A::Matrix{Float64}, b::Vector{Float64}
	)
		n = size(A, 1)
		L = zeros(Float64, n, n)

		for i in 1:n
			sum_square = sum(L[i, 1:i-1] .^ 2)
			L[i, i] = sqrt(A[i, i] - sum_square)

			for j in i+1:n
				sum_ele = sum(L[i, 1:i-1] .* L[j, 1:i-1])
				L[j, i] = (A[i, j] - sum_ele) / L[i, i]
			end
		end

		y = zeros(Float64, n)
		for i in 1:n
		    y[i] = (b[i] - dot(L[i,1:i-1], y[1:i-1])) / L[i,i]
		end

		x = zeros(Float64, n)
		for i in n:-1:1
		    x[i] = (y[i] - dot(L[1:i-1,i], x[1:i-1])) / L[i,i]
		end

		println(x)
	end

	B = [4.0 2.0 2.0;
	     2.0 10.0 5.0;
	     2.0 5.0 9.0]
	
	cholesky(B, b)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.2"
manifest_format = "2.0"
project_hash = "f352ceee806168c8ae38887a01d7bae6ca62470b"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╟─6ca1508e-d9b3-11f0-85e3-a9cca3c67ff7
# ╠═61958f32-a8b6-40d6-be2e-858574ab7708
# ╟─f880da5a-0dd1-4592-bf18-3d830441abd4
# ╠═41c19325-2ac9-42f7-862e-37abd505bfc8
# ╟─ce647298-2958-48fa-9361-3f48955a098f
# ╠═317227db-1cf5-4329-8441-9ad00a9b1013
# ╟─85485111-3749-4459-9b8b-d6b3186435e5
# ╠═16e9814d-42df-4fa8-8b44-e81565e3d335
# ╟─7aea5bfd-3757-48fe-adc7-4860e3bc68a8
# ╠═beddc0f4-18c0-4b22-aad8-d173023dc259
# ╟─cb4087f2-5021-4f48-b244-1475303bfd2d
# ╠═2e7b26e0-c6be-46f0-bcd2-2cd092620217
# ╟─09cf821b-991e-4e5c-823b-726a25414d9c
# ╠═f97c4dab-dba6-4038-bd18-64c959a2a8e2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
