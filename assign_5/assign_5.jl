using Plots
using HermiteInterpolation
using DataInterpolations

# Task 1
function divided_differences(x, y)
    n = length(x)
    a = copy(y)
    for j in 2:n
        for i in n:-1:j
            a[i] = (a[i] - a[i-1]) / (x[i] - x[i-j+1])
        end
    end
    return a
end

function newton_evaluate(x_data, a, x)
    n = length(a)
    result = a[n]
    for i in n-1:-1:1
        result = result * (x - x_data[i]) + a[i]
    end
    return result
end

function find_root(f, a, b, target; tol=1e-6)
    while abs(b - a) > tol
        mid = (a + b) / 2
        if f(mid) < target
            a = mid
        else
            b = mid
        end
    end
    return (a + b) / 2
end

x_vals = [-1.0, 0.0, 2.0, 3.0]
y_vals = [-4.0, -1.0, 0.0, 3.0]

a = divided_differences(x_vals, y_vals)
f_interp = x -> newton_evaluate(x_vals, a, x)

println("f(1.5) ≈ ", f_interp(1.5))

x_root = find_root(f_interp, 0.0, 3.0, 0.5)
println("f(x) = 0.5 at x ≈ ", x_root)

x_vals_new = [-1.0, 0.0, 2.0, 3.0, 3.1]
y_vals_new = [-4.0, -1.0, 0.0, 3.0, 0.1]
a_new = divided_differences(x_vals_new, y_vals_new)
f_interp_new = x -> newton_evaluate(x_vals_new, a_new, x)

xs = -1.5:0.01:3.5
ys_old = [f_interp(x) for x in xs]
ys_new = [f_interp_new(x) for x in xs]

p1 = plot(xs, ys_old, label="Newton (4 points)", lw=2, color=:red)
plot!(p1, xs, ys_new, label="Newton (5 points)", lw=2, ls=:dash)
scatter!(p1, x_vals, y_vals, label="Original Data", color=:black)
scatter!(p1, [3.1], [0.1], label="Added Point", color=:red)

# Other methods:

# Hermite Interpolation:
derivative_set = [3.0, 0.5, 3.0, -1.0]
hermite = HermiteInterpolation.fit(x_vals, y_vals, derivative_set)
x_hermite = -1.5:0.01:3.5
y_hermite = hermite.(x_hermite)
plot!(p1, x_hermite, y_hermite, label="Hermite", lw=2, color=:green)

# LaGrange Interpolation:
lagrange = LagrangeInterpolation(
    y_vals, x_vals; extrapolation = ExtrapolationType.Linear
)

x_lagrange = -1.5:0.01:3.5
y_lagrange = lagrange.(x_lagrange)
p2 = plot(
    x_lagrange, y_lagrange, label="LaGrange",
    lw=2, color=:blue, ls=:dashdot
)

# Cubic Spline Interpolation:
cubic_spline = CubicSpline(
    y_vals, x_vals; extrapolation = ExtrapolationType.Linear
)

x_cspline = -1.5:0.01:3.5
y_cspline = lagrange.(x_lagrange)
plot!(
    p2, x_cspline, y_cspline, label="Cubic-Spline",
    lw=2, color=:cyan, ls=:dot
)

# Task 2
x_2 = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
y_2 = [1/26, 1/17, 1/10, 1/5, 1/2, 1.0, 1/2, 1/5, 1/10, 1/17]

a_2 = divided_differences(x_2, y_2)
f_interp_2 = x -> newton_evaluate(x_2, a_2, x)
xs_2 = -5.5:0.01:4.5
ys_2 = [f_interp_2(x) for x in xs_2]

p3 = plot(xs_2, ys_2, label="Newton Task 2", lw=2, color=:red)
scatter!(p3, x_2, y_2, label="Original Data", color=:black)

plot(p1, p2, layout=(2,1))
