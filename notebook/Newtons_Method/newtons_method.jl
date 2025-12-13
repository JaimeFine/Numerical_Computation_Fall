using Printf

# Define the function and its derivative
f(x) = x^2 - 2      # An example: root of x^2 is √2
df(x) = 2x          # Its derivative

# Newton-Raphson method:
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

root = newton(f, df, 1.0)
@printf("Estimated root: %.8f\n", root)

# For default f64 precision (16-digits):
# println("Estimated root: ", root)