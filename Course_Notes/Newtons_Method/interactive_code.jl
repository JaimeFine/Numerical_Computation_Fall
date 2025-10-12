using Printf
# Note that this code is different than the newtons_method.jl
# You can apply either, both are correct!
println("Enter the function:")
user_input_f = readline()
f = eval(Meta.parse("x -> $user_input_f"))

println("Enter its derivative:")
user_input_df = readline()
df = eval(Meta.parse("x -> $user_input_df"))

println("Enter the initial value:")
x0 = parse(Float64, readline())

println("Enter your precision:")
e = parse(Float64, readline())

println("Enter the max iteration:")
max_iter = parse(Int, readline())

function newton(f, df, x0, e, max_iter)
    flag = false
    count = 0
    x = x0
    x_new = x0

    while true
        fx = f(x)
        dfx = df(x)

        if abs(dfx) < eps()
            error("Derivative too small - possible division by zero.")
        end

        x_new = x - fx / dfx
        count += 1

        if abs(x_new - x) <= e
            flag = true
            break
        elseif count <= max_iter
            x = x_new
        else
            break
        end
    end

    if flag
        d = max(0, -floor(Int, log10(e)))
        format = Printf.Format("The root is: %." * string(d) * "f\n")
        result = Printf.format(format, x_new)
        print(result)
    else
        println("Failed to converge within $max_iter iterations...")
    end
end

newton(f, df, x0, e, max_iter)