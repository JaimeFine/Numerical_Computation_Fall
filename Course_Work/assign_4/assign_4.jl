using Printf
using BenchmarkTools

f(x) = cos(x) + 1.0 / (1.0 + exp(-2.0 * x))
df(x) = -sin(x) + 2.0 * exp(-2.0 * x) / (1.0 + exp(-2.0 * x))^2

function iteration_1(x0; max_iter=100, tol=1e-8)
    x = x0
    for i in 1:max_iter
        x_new = acos(-1.0 / (1.0 + exp(-2 * x)))

        if abs(x_new - x) < tol
            @printf(
                "The minimum positive root of the first scheme is: %.5f",
                x_new
            )
            return x_new
        end
        x = x_new
    end
    println("Unable to solve after $max_iter iterations.")
    return x_new
end

function iteration_2(x0; max_iter=100, tol=1e-8)
    x = x0
    for i in 1:max_iter
        arg = -1.0 / (1.0 + 1.0 / cos(x))
        if arg <= 0.0
            @printf("Error! Log of non-positive number at x = %.5f.", x)
            return nothing
        end

        x_new = 0.5 * log(arg)

        if abs(x_new - x) < tol
            println(
                "The minimum positive root of the second scheme is: %.5f",
                x_new
            )
            return x_new
        end
        x = x_new
    end
    println("Unable to solve after $max_iter iterations.")
    return x_new
end

function newton(f, df, x0; max_iter=100, tol=1e-8)
    x = x0
    for i in 1:max_iter
        fx = f(x)
        dfx = df(x)
        if abs(dfx) < eps()
            error("Derivative too small...")
        end

        x_new = x - fx / dfx

        if abs(x_new - x) < tol
            @printf(
                "The minimum positive root of the newton's method is: %.5f",
                x_new
            )
            return x_new
        end
        x = x_new
    end
    println("Unable to solve after $max_iter iterations.")
    return x_new
end

@time iteration_1(3.0)
println("-----------------------------------")
iteration_2(3)
println("\n-----------------------------------")
@time newton(f, df, 3.0)