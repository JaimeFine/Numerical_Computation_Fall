using LinearAlgebra
using Statistics

function estimate_rigid_transform(P::AbstractMatrix{<:Real}, Q::AbstractMatrix{<:Real})
    @assert size(P,1) == 3 && size(Q,1) == 3 && size(P,2) == size(Q,2)

    μP = mean(P, dims=2)
    μQ = mean(Q, dims=2)

    Pc = P .- μP
    Qc = Q .- μQ

    H = Pc * Qc'
    s = svd(H)
    U = s.U
    S = s.S
    V = s.V

    R = V * U'

    if det(R) < 0
        # flip the last column of V (handle general size with `end`)
        V[:, end] .= -V[:, end]
        R = V * U'
    end

    t = μQ .- R * μP
    return R, t
end

function rodrigues_from_R(R::AbstractMatrix{<:Real})
    θ = acos(clamp((tr(R) - 1)/2, -1, 1))
    if θ < 1e-12
        return zeros(3), 0.0
    end
    v = [R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]] / (2*sin(θ))
    r = v * θ
    return r, θ
end


# Initial 3d dots (homogeneous rows):
P_raw = []

# Dots after the transformation (homogeneous rows):
Q_raw = []

# Convert the provided row-form homogeneous data to 3×n point matrices
# (columns are points) by taking the first 3 columns and transposing.
P = permutedims(P_raw[:, 1:3])  # 3×n
Q = permutedims(Q_raw[:, 1:3])  # 3×n

R, t = estimate_rigid_transform(P, Q)
r, θ = rodrigues_from_R(R)

# Expected answers:
# Rotation: a, b, c = 25.399051233106778 -54.101700051821325   5.637043873925570
# Panning: x, y, z = -15.174821832828904   9.590334711153755  21.176335086376653

rot = r * 180 / pi
println("Rodriguez's vector r: ", rot)
println("Translation t: ", vec(t))
