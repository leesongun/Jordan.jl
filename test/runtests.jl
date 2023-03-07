using Jordan
using Test
using Symbolics
using LinearAlgebra

Oct(x, y) = all(isequal(simplify((x-y).v[i]; expand=true), 0) for i in 1:8)
makeOct(x) = Octonion([Symbolics.scalarize(x[i]) for i in 1:8])

@testset "Octonions" begin
    @variables x[1:3, 1:8]
    a, b, c = [makeOct(x[i, :]) for i in 1:3]

    # addition: abelian group
    # unital magma
    @test Oct(a + 0, a)
    @test Oct(0 + a, a)

    @test Oct(-a + a, 0)
    @test Oct(a + -a, 0)

    @test Oct((a + b) + c, a + (b + c))
    @test Oct(a + b, b + a)

    # multiplication: Moufang loop
    @test Oct(a * 1, a)
    @test Oct(1 * a, a)

    @test Oct(a * (inv(a) * b), b)
    @test Oct((b * inv(a)) * a, b)

    @test Oct((a * b * a) * c, a * (b * (a * c))) # left Bol
    @test Oct(c * (a * b * a), ((c * a) * b) * a) # right Bol, redundant due to conj
    @test Oct(a * (b * c) * a, (a * b) * (c * a)) # Moufang identity

    # conjugation: involutive anti-automorphism
    @test Oct(a, conj(conj(a)))

    @test Oct(conj(Octonion(0)), 0)
    @test Oct(conj(Octonion(1)), 1)

    @test Oct(conj(a) + conj(b), conj(b + a))
    @test Oct(conj(a) * conj(b), conj(b * a))

    # algebra: distributive
    @test Oct((a + b) * c, a * c + b * c)
    @test Oct(a * (b + c), a * b + a * c) # redundant due to conj anti-automorphism
end

Alb(x, y) = all(Oct.(x - y, 0))
makeAlb(x) = Hermitian(reshape([makeOct(x[i, :]) for i in 1:9], (3, 3)))

@testset "Albert" begin
    @variables x[1:2, 1:9, 1:8]

    a, b = [makeAlb(x[i, :, :]) for i in 1:2]
    t(x, y) = (x * y + y * x)# / 2
    @test Alb(t(a, b), t(b, a))
    @test Alb(t(t(a, b), t(a, a)), t(a, t(b, t(a, a))))
end
