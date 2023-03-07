using Jordan
using Test
using Symbolics

function Oct(x, y)
    z = (x - y).v
    return all(isequal(simplify(z[i]; expand=true), 0) for i in 1:8)
end

@testset "Octonions" begin
    @variables x[1:3, 1:8]
    (a, b, c) = [Octonion([Symbolics.scalarize(x[i, j]) for j in 1:8]) for i in 1:3]

    @test Oct(a, conj(conj(a)))

    @test Oct(a + b, b + a)
    @test Oct(conj(a) * conj(b), conj(b * a))

    @test Oct((a + b) * c, a * c + b * c)
    @test Oct(a * (b + c), a * b + a * c) # redundant due to conj anti-automorphism

    # Alternative algebra (any 2 generates associative algebra)
    @test Oct((a * a) * b, a * (a * b))
    @test Oct((a * b) * b, a * (b * b))
    @test Oct((a * b) * a, a * (b * a)) # redundant due to distribution law

    # Moufang loop
    @test Oct((a * (b * a)) * c, a * (b * (a * c))) # left Bol
    @test Oct(c * ((a * b) * a), ((c * a) * b) * a) # right Bol, redundant due to conj anti-automorphism
    @test Oct((a * b) * (c * a), (a * (b * c)) * a) # Moufang identity. redundant by definition.
    @test Oct((a * b) * (c * a), a * ((b * c) * a)) # Moufang identity. redundant by definition and flexibility.
end

function MakeAlbert(x, i, j)
    if i < j
        return Octonion([Symbolics.scalarize(x[6-i-j, k]) for k in 1:8])
    elseif i == j
        return Octonion(x[i, 9])
    else
        return conj(MakeAlbert(x, j, i))
    end
end

function MakeMatrix(x)::Matrix
    return reshape([MakeAlbert(x, i, j) for (i, j) in Iterators.product(1:3, 1:3)], (3, 3))
end

function Alb(x, y)
    t(a) = Oct(a, 0)
    return all(t.(x - y))
end

@testset "Albert" begin
    @variables x[1:2, 1:3, 1:9]

    (a, b) = [MakeMatrix(x[i, :, :]) for i in 1:2]
    t(x, y) = (x * y + y * x)# / 2
    @test Alb(t(a, b), t(b, a))
    @test Alb(t(t(a, b), t(a, a)), t(a, t(b, t(a, a))))
end
