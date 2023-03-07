module Jordan
using StaticArrays

export Octonion, Jordan

struct Octonion{T<:Real} <: Number
    v::SVector{8,T}
end

Octonion(s::Vector{Real}) = Octonion(SVector{8}(s))
Octonion{T}(x::Real) where {T<:Real} = Octonion(convert(T, x))
Octonion{T}(o::Octonion) where {T<:Real} = Octonion{T}(o.v)
Octonion(x::Real) = Octonion(hcat(x, zeros(typeof(x), 7)))

Base.promote_rule(::Type{Octonion{T}}, ::Type{S}) where {T<:Real,S<:Real} = Octonion{promote_type(T, S)}
Base.promote_rule(::Type{Octonion{T}}, ::Type{Octonion{S}}) where {T<:Real,S<:Real} = Octonion{promote_type(T, S)}

Base.:-(o::Octonion) = Octonion(-o.v)

Base.:+(o1::Octonion, o2::Octonion) = Octonion(o1.v + o2.v)
Base.:-(o1::Octonion, o2::Octonion) = Octonion(o1.v - o2.v)
Base.conj(o::Octonion) = Octonion(vcat(o.v[1], -o.v[2:8]))

# Warning: not numerically best way to implement!
function Base.:*(o1::Octonion, o2::Octonion)::Octonion
    i1::SVector{7} = o1.v[2:8]
    i2::SVector{7} = o2.v[2:8]
    rotmul(s) = circshift(i1, s) .* i2 - i1 .* circshift(i2, s)
    mul = sum(circshift(rotmul(a), 2 * a) for a in [1, 2, 4])
    ret = o1.v[1] * o2.v + vcat(-i1'i2, i1 * o2.v[1] + mul)
    return Octonion(ret)
end

Base.:(==)(o1::Octonion, o2::Octonion) = o1.v == o2.v
Base.isequal(o1::Octonion, o2::Octonion) = isequal(o1.v, o2.v)

struct Albert{T}
    m::SMatrix{Octonion{T},3,3,9}
    # inner constructor to enforce adjoint condition
end
Base.:*(a1::Albert, a2::Albert)::Albert = (a1.m * a2.m + a2.m * a1.m) / 2

end
