using Plots

x = [1.0, 2.0, 3.0, 4.0]
y = [60.0, 30.0, 20.0, 15.0]
# Transform it into linear equation:
# ln(y) = bx + ln(a)
y_new = log.(y)

function build_matrix(x)
    n = length(x)
    X = Array{Float64}(undef, n, 2)
    for i in 1:n
        X[i, 1] = 1.0
        X[i, 2] = x[i]
    end
    return X
end

X = build_matrix(x)

Xt = transpose(X)
XtX = Xt * X
Xty = Xt * y_new

XtX_inv = inv(XtX)
result = XtX_inv * Xt * y_new
b = result[2]
a = exp(result[1])

println("a is ", a, ", b is ", b)

x_r = -10.0:0.01:20.0
y_r = [a * exp(b * xi) for xi in x_r]

plot(x_r, y_r)
