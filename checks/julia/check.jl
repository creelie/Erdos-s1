# check.jl -- a third, independent check in Julia with BigInt arithmetic.
#
#   julia checks/julia/check.jl
#
# For both sets of the paper it checks general position, every relation of the
# certificate (positive weights, the identity, the determinant and the Cramer numerators),
# the covering condition, and, by a direct exact test that ignores the certificate, that
# no 8 points are in convex position.  Only the Julia standard library is used; the JSON
# files are read as flat lists of integers.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

ints(path) = [parse(BigInt, m.match) for m in eachmatch(r"-?\d+", read(path, String))]

function bareiss(A::Matrix{BigInt})
    A = copy(A); n = size(A, 1); s = 1; prev = big(1)
    for i in 1:n-1
        if A[i, i] == 0
            r = findfirst(r -> A[r, i] != 0, i+1:n)
            r === nothing && return big(0)
            r += i
            A[[i, r], :] = A[[r, i], :]; s = -s
        end
        for r in i+1:n, c in i+1:n
            A[r, c] = div(A[r, c] * A[i, i] - A[r, i] * A[i, c], prev)
        end
        prev = A[i, i]
    end
    return s * A[n, n]
end

orient(P, idx) = bareiss(BigInt[j == 1 ? 1 : P[i][j-1] for i in idx, j in 1:length(P[1])+1])

function combos(n, k)
    out = Vector{Vector{Int}}()
    c = collect(1:k)
    while true
        push!(out, copy(c))
        i = k
        while i >= 1 && c[i] == n - k + i
            i -= 1
        end
        i == 0 && return out
        c[i] += 1
        for j in i+1:k
            c[j] = c[j-1] + 1
        end
    end
end

function inside(P, p, S)
    o = orient(P, S)
    o == 0 && return false
    for j in eachindex(S)
        T = copy(S); T[j] = p
        t = orient(P, T)
        (t == 0 || (t > 0) != (o > 0)) && return false
    end
    return true
end

function check(name, wfile, cfile, d, k, nrel)
    v = ints(joinpath(ROOT, "highdim", wfile))
    P = [v[i*d+1:i*d+d] for i in 0:div(length(v), d)-1]
    n = length(P)
    c = ints(joinpath(ROOT, "highdim", cfile))
    rec = 2d + 4
    @assert length(c) % rec == 0
    passed = 0
    gp = all(orient(P, S) != 0 for S in combos(n, d + 1))
    @assert gp "$name: general position"
    passed += 1
    supports = Vector{Set{Int}}()
    for r in 0:div(length(c), rec)-1
        f = c[r*rec+1:(r+1)*rec]
        p = Int(f[1]) + 1
        S = [Int(x) + 1 for x in f[2:d+2]]
        m = f[d+3:2d+3]
        M = sum(m)
        @assert all(m .> 0) "$name: weights"
        @assert all(sum(m[j] * P[S[j]][t] for j in 1:d+1) == M * P[p][t] for t in 1:d) "$name: identity"
        D = orient(P, S)
        @assert D == f[end] && D != 0 && mod(D, M) == 0 "$name: determinant"
        for j in 1:d+1
            T = copy(S); T[j] = p
            N = orient(P, T)
            @assert m[j] * D == M * N "$name: Cramer numerator"
        end
        push!(supports, Set([p; S]))
        passed += 1
    end
    @assert length(supports) == nrel
    for T in combos(n, k)
        Ts = Set(T)
        @assert any(issubset(s, Ts) for s in supports) "$name: uncovered subset $T"
    end
    passed += 1
    gons(m) = count(T -> !any(inside(P, T[a], S) for a in eachindex(T)
                              for S in (T[setdiff(eachindex(T), a)][s] for s in combos(m - 1, d + 1))),
                    combos(n, m))
    g8 = gons(k)
    @assert g8 == 0 "$name: $g8 convex $k-subsets"
    g7 = gons(k - 1)
    @assert g7 > 0
    passed += 1
    println("$name: $n points in R^$d, $(length(supports)) relations valid, every $k-subset covered, ",
            "convex $k-subsets: $g8 (control: convex $(k-1)-subsets: $g7)")
    return passed
end

total = check("X4", "witness_r4_11.json", "cert_r4_11.json", 4, 8, 28) +
        check("X3", "witness_r3_15.json", "cert_r3_15.json", 3, 8, 161)
println("$total checks passed, 0 failed")
