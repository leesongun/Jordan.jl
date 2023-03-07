module Jordan
using StaticArrays

export Octonion, Jordan

struct Octonion{T<:Real} <: Number
    v::SizedVector{8,T}
end

# temporary solution
Octonion(s::Vector) = Octonion(SizedVector{8}(s))
Octonion(s::Vector{Real}) = Octonion(SizedVector{8}(s))
Octonion{T}(x::Real) where {T<:Real} = Octonion(convert(T, x))
Octonion{T}(o::Octonion) where {T<:Real} = Octonion{T}(o.v)
Octonion(x::Real) = Octonion(vcat(x, zeros(typeof(x), 7)))

Base.promote_rule(::Type{Octonion{T}}, ::Type{S}) where {T<:Real,S<:Real} = Octonion{promote_type(T, S)}
Base.promote_rule(::Type{Octonion{T}}, ::Type{Octonion{S}}) where {T<:Real,S<:Real} = Octonion{promote_type(T, S)}

Base.:-(o::Octonion) = Octonion(-o.v)
Base.real(o::Octonion) = o.v[1]
Base.conj(o::Octonion) = Octonion(vcat(o.v[1], -o.v[2:8]))
Base.abs2(o::Octonion) = o.v'o.v
Base.inv(o::Octonion) = Octonion(conj(o).v ./ abs2(o))

Base.:+(o1::Octonion, o2::Octonion) = Octonion(o1.v + o2.v)
Base.:-(o1::Octonion, o2::Octonion) = Octonion(o1.v - o2.v)

# Warning: not numerically best way to implement!
function Base.:*(o1::Octonion, o2::Octonion)::Octonion
    i1, i2 = o1.v[2:8], o2.v[2:8]
    rotmul(s) = circshift(circshift(i1, s) .* i2 - i1 .* circshift(i2, s), 2 * s)
    return Octonion(o1.v[1] * o2.v + vcat(-i1'i2, i1 * o2.v[1] + sum(rotmul.([1, 2, 4]))))
end

Base.:(==)(o1::Octonion, o2::Octonion) = o1.v == o2.v
Base.isequal(o1::Octonion, o2::Octonion) = isequal(o1.v, o2.v)

end
