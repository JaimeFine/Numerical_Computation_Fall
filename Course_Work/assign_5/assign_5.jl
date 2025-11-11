using Plots

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

plot(xs, ys_old, label="Newton (4 points)", lw=2)
plot!(xs, ys_new, label="Newton (5 points)", lw=2, ls=:dash)
scatter!(x_vals, y_vals, label="Original Data", color=:black)
scatter!([3.1], [0.1], label="Added Point", color=:red)