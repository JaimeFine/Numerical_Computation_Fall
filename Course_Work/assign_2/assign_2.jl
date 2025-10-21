using LinearAlgebra

function jacobi(A::Matrix{Float64}, b::Vector{Float64}; x0=zeros(size(b)), tol=1e-8, max_iter=1000)
    n = length(b)
    x = copy(x0)
    D = Diagonal(diag(A))
    L = LowerTriangular(A) - D
    U = UpperTriangular(A) - D

    D_inv = inv(D)

    for k in 1:max_iter
        x_new = -D_inv * (L + U) * x + D_inv * b
        if norm(x_new - x, Inf) < tol
            return x_new
        end
        x = x_new
    end

    return x
end

function gauss_seidel(A::Matrix{Float64}, b::Vector{Float64}; x0=zeros(size(b)), tol=1e-8, max_iter=1000)
    n = length(b)
    x = copy(x0)

    for k in 1:max_iter
        for i in 1:n
            sum1 = sum((A[i, j] * x[j] for j in 1:i-1), init=0.0)
            sum2 = sum((A[i, j] * x[j] for j in i+1:n), init=0.0)
            x[i] = (b[i] - sum1 - sum2) / A[i, i]
        end
        if norm(A * x - b, Inf) < tol
            return x
        end
    end

    return x
end

function sor(A::Matrix{Float64}, b::Vector{Float64}, w::Float64; x0=zeros(size(b)), tol=1e-8, max_iter=1000)
    n = length(b)
    x = copy(x0)

    for k in 1:max_iter
        for i in 1:n
            sum1 = sum((A[i, j] * x[j] for j in 1:i-1), init=0.0)
            sum2 = sum((A[i, j] * x[j] for j in i+1:n), init=0.0)
            x[i] = (1 - w) * x[i] + (w / A[i, i]) * (b[i] - sum1 - sum2)
        end
        if norm(A * x - b, Inf) < tol
            return x
        end
    end

    return x
end


A = [4.0 3.0 0.0;
     3.0 4.0 -1.0;
     0.0 -1.0 4.0]

b = [24.0, 30.0, -24.0]

x = jacobi(A, b)
println("Solution Jacobi: ", x)

x = gauss_seidel(A, b)
println("Solution Gauss-Seidel: ", x)

x = sor(A, b, 1.8)
println("Solution SOR with w = 1.8: ", x)

x = sor(A, b, 1.22)
println("Solution SOR with w = 1.22: ", x)