# Note that this code is different than the newtons_method.jl
# You can apply either, both are correct!
println("Enter the function:")
f = eval(Meta.parse("x -> $readline()"))
println("Enter its derivative:")
df = eval(Meta.parse("x -> $readline()"))
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
    adj_e = (1 - abs(df(x0))) * e

    while true
        fx = f(x)
        dfx = df(x)

        if abs(dfx) < eps()
            error("Derivative too small - possible division by zero.")
        end

        x_new = x - fx / dfx
        count += 1

        if abs(x_new - x) <= adj_e
            flag = true
            break
        elseif count <= max_iter
            x = x_new
        else
            break
        end
    end

    if flag
        println("The root is: ", x_new)
    else
        println("Failed to converge within $max_iter iterations...")
    end
end

newton(f, df, x0, e, max_iter)