using Printf

# Newton-Raphson method
function newton(f, df, x0; e=1e-8, max_iter=100)
    x = x0
    for i in 1:max_iter
        fx = f(x)
        dfx = df(x)
        if abs(dfx) < eps()
            error("Derivative too small - possible division by zero.")
        end
        x_new = x - fx / dfx
        if abs(x_new - x) < e
            return x_new
        end
        x = x_new
    end
    error("Did not succeed to converge within $max_iter iterations.")
end

# Problem 1
f1(x) = x^3 - x^2 - x - 1
df1(x) = 3x^2 - 2x - 1
root1 = newton(f1, df1, 2.0)
@printf("Root of x^3 - x^2 - x - 1 = 0: %.8f\n", root1)

# Problem 2
f2(x) = cos(x) - 0.5 - sin(x)
df2(x) = -sin(x) - cos(x)
root2 = newton(f2, df2, 1.0)
@printf("Smallest positive root of cos(x) = 0.5 + sin(x): %.8f\n", root2)
