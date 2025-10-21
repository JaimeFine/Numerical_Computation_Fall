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

A = [10.0 -1.0 2.0 0.0;
     -1.0 11.0 -1.0 3.0;
     2.0 -1.0 10.0 -1.0;
     0.0 3.0 -1.0 8.0]

b = [6.0, 25.0, -11.0, 15.0]

x = jacobi(A, b)
println("Solution: ", x)
