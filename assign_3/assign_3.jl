using LinearAlgebra
using BenchmarkTools
using Random

# Generate the data:
const N = 5000
Random.seed!(42)

A = rand(N, N)
x_true = rand(N)
b = A * x_true

# 1: Direct library solve (Time = T1):
print("[1: Direct library solve]")
@time x_lib = A \ b

T1 = @benchmark $A \ $b samples=2 seconds=30

# 2: LU Decomposition (Time T2):
function lu_solve(A, b)
    F_lu = lu(A)
    x = F_lu \ b
    return x
end

print("[2: LU decomposition solve]")
@time x_LU = lu_solve(A, b)

T2 = @benchmark lu_solve($A, $b) samples=2 seconds=30

# norm(A*X-b) verification:
residual_norm_lib = norm(A * x_lib - b)
residual_norm_LU = norm(A * x_LU - b)
initial_norm = norm(A * x_true - b)

println("Verification:")
println("Initial norm (error check on b): $(initial_norm)")
println("Residual norm (Direct Solve):    $(residual_norm_lib)")
println("Residual norm (LU Solve):        $(residual_norm_LU)")

println("================================================")

println(">>> Comparison:")
println("T1 (Direct Solve):")
display(T1)
println("T2 (LU Decomposition + Solve):")
display(T2)

println("================================================")

println("Second exercise:")
A_spec = [
    0.8147   0.0975   0.1576   0.1419   0.6557;
    0.9058   0.2785   0.9706   0.4218   0.0357;
    1.27e9   0.5469   0.9572   0.9157   0.8491;
    0.9134   0.9575   4.854e7  0.7922   0.9340;
    0.6324   0.9649   0.8003   0.9595   0.6787
]

b_spec = [
    2.258;
    1.5977;
    1.2700000023549e9;
    2.42700039042e7;
    3.36025
]

x_spec = A_spec \ b_spec

# Verification:
residual_norm_spec = norm(A_spec * x_spec - b_spec)

println("Solution vector x:")
display(x_spec)

println("Verification Residual Norm ||Ax - b||:")
println("$(residual_norm_spec)")