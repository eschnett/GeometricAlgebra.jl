using GeometricAlgebra
using Base.Test

using FlexibleArrays
using Quaternions

typealias Real0f MultiVector(0, Float64)
typealias Complex1f MultiVector(1, Float64)
typealias Quaternion2f MultiVector(2, Float64)
typealias MultiVector2i MultiVector(2, Int)
typealias MultiVector3f MultiVector(3, Float64)
typealias MultiVector3i MultiVector(3, Int)

typealias MultiVector4f let
    D = 4
    Cycles = (2,2,2,2)
    Parities = (+1,-1,-1,-1)
    Commutators = ((0,-1,-1,-1), (-1,0,-1,-1), (-1,-1,0,-1), (-1,-1,-1,-1))
    MultiVector(D, Cycles, Parities, Commutators, Float64)
end

# typealias MultiVector9f MultiVector(9, Float64)

@test length(Real0f) == 2^0
@test length(Complex1f) == 2^1
@test length(MultiVector2i) == 2^2
@test length(MultiVector3f) == 2^3
@test length(MultiVector4f) == 2^4
# @test length(MultiVector9f) == 2^9

m0 = Real0f()
m1 = Complex1f()
m2 = MultiVector2i()
m3 = MultiVector3f()
m4 = MultiVector4f()
# m9 = MultiVector9f()

@test length(m0) == 2^0
@test length(m1) == 2^1
@test length(m2) == 2^2
@test length(m3) == 2^3
@test length(m4) == 2^4
# @test length(m9) == 2^9

@test length(collect(m0)) == length(m0)
@test length(collect(m1)) == length(m1)
@test length(collect(m2)) == length(m2)
@test length(collect(m3)) == length(m3)
@test length(collect(m4)) == length(m4)
# @test length(collect(m9)) == length(m9)

@test typeof(collect(m0)[1]) == eltype(m0)
@test typeof(collect(m1)[1]) == eltype(m1)
@test typeof(collect(m2)[1]) == eltype(m2)
@test typeof(collect(m3)[1]) == eltype(m3)
@test typeof(collect(m4)[1]) == eltype(m4)
# @test typeof(collect(m9)[1]) == eltype(m9)

@test all(x->x==0, m0)
@test all(x->x==0, m1)
@test all(x->x==0, m2)
@test all(x->x==0, m3)
@test all(x->x==0, m4)
# @test all(x->x==0, m9)

m0 = setindex(m0, 1)

m1 = setindex(m1, 1, 0)
m1 = setindex(m1, 2, 1)

m2 = setindex(m2, 1, 0,0)
m2 = setindex(m2, 2, 0,1)
m2 = setindex(m2, 3, 1,0)
m2 = setindex(m2, 4, 1,1)

m3 = setindex(m3, 1, 0,0,0)
m3 = setindex(m3, 2, 0,0,1)
m3 = setindex(m3, 3, 0,1,0)
m3 = setindex(m3, 4, 0,1,1)
m3 = setindex(m3, 5, 1,0,0)
m3 = setindex(m3, 6, 1,0,1)
m3 = setindex(m3, 7, 1,1,0)
m3 = setindex(m3, 8, 1,1,1)

let
    count = 0
    for i in CartesianRange(CartesianIndex(0,0,0,0), CartesianIndex(1,1,1,1))
        count += 1
        m4 = setindex(m4, count, i)
    end
end

# let
#     count = 0
#     for i in CartesianRange(CartesianIndex(0,0,0,0,0,0,0,0,0),
#                             CartesianIndex(1,1,1,1,1,1,1,1,1))
#         count += 1
#         m9 = setindex(m9, count, i)
#     end
# end

@test m0[] == 1

@test m1[0] == 1
@test m1[1] == 2

@test m2[0,0] == 1
@test m2[0,1] == 2
@test m2[1,0] == 3
@test m2[1,1] == 4

@test m3[0,0,0] == 1
@test m3[0,0,1] == 2
@test m3[0,1,0] == 3
@test m3[0,1,1] == 4
@test m3[1,0,0] == 5
@test m3[1,0,1] == 6
@test m3[1,1,0] == 7
@test m3[1,1,1] == 8

let
    count = 0
    for i in CartesianRange(CartesianIndex(0,0,0,0), CartesianIndex(1,1,1,1))
        count += 1
        @test m4[i] == count
    end
end

# let
#     count = 0
#     for i in CartesianRange(CartesianIndex(0,0,0,0,0,0,0,0,0),
#                             CartesianIndex(1,1,1,1,1,1,1,1,1))
#         count += 1
#         @test m9[i] == count
#     end
# end

z0 = zero(Real0f)
z1 = zero(Complex1f)
z2 = zero(MultiVector2i)
z3 = zero(MultiVector3f)
z4 = zero(MultiVector4f)
# z9 = zero(MultiVector9f)

@test z0[] == 0

@test z1[0] == 0
@test z1[1] == 0

@test z2[0,0] == 0
@test z2[0,1] == 0
@test z2[1,0] == 0
@test z2[1,1] == 0

for i in CartesianRange(CartesianIndex(0,0,0), CartesianIndex(1,1,1))
    @test z3[i] == 0
end

for i in CartesianRange(CartesianIndex(0,0,0,0), CartesianIndex(1,1,1,1))
    @test z4[i] == 0
end

# for i in CartesianRange(CartesianIndex(0,0,0,0,0,0,0,0,0),
#                         CartesianIndex(1,1,1,1,1,1,1,1,1))
#     @test z9[i] == 0
# end



s0 = 2.0 * m0
s1 = 2.0 * m1
s2 = 2.0 * m2
s3 = 2.0 * m3
s4 = 2.0 * m4
# s9 = 2.0 * m9

s0 = +s0 + m0
s1 = +s1 + m1
s2 = +s2 + m2
s3 = +s3 + m3
s4 = +s4 + m4
# s9 = +s9 + m9

s0 = s0 - -m0
s1 = s1 - -m1
s2 = s2 - -m2
s3 = s3 - -m3
s4 = s4 - -m4
# s9 = s9 - -m9

s0 = s0 * (1/4)
s1 = s1 * (1/4)
s2 = s2 * (1/4)
s3 = s3 * (1/4)
s4 = s4 * (1/4)
# s9 = s9 * (1/4)

@test s0[] == 1

@test s1[0] == 1
@test s1[1] == 2

@test s2[0,0] == 1
@test s2[0,1] == 2
@test s2[1,0] == 3
@test s2[1,1] == 4

@test s3[0,0,0] == 1
@test s3[0,0,1] == 2
@test s3[0,1,0] == 3
@test s3[0,1,1] == 4
@test s3[1,0,0] == 5
@test s3[1,0,1] == 6
@test s3[1,1,0] == 7
@test s3[1,1,1] == 8

let
    count = 0
    for i in CartesianRange(CartesianIndex(0,0,0,0), CartesianIndex(1,1,1,1))
        count += 1
        @test s4[i] == count
    end
end

# let
#     count = 0
#     for i in CartesianRange(CartesianIndex(0,0,0,0,0,0,0,0,0),
#                             CartesianIndex(1,1,1,1,1,1,1,1,1))
#         count += 1
#         @test s9[i] == count
#     end
# end

@test Real0f(2.0) * Real0f(3.0) == Real0f(2.0 * 3.0)

let
    c = Complex(2.0 + 3.0im) * Complex(4.0 + 5.0im)
    @test (Complex1f(2.0, 3.0) * Complex1f(4.0, 5.0) ==
           Complex1f(real(c), imag(c)))
end

let
    q = Quaternion(2.0,3.0,4.0,5.0) * Quaternion(6.0,7.0,8.0,9.0)
    @test (Quaternion2f(2.0,3.0,4.0,5.0) * Quaternion2f(6.0,7.0,8.0,9.0) ==
           Quaternion2f(q.s, q.v1, q.v2, q.v3))
end


let
    v2m(vals) = MultiVector3f(vals[1],0.0,0.0,vals[2], 0.0,vals[3],vals[4],0.0)
    q2v(q) = (qr.s, qr.v1,qr.v2,qr.v3)
    q1 = Quaternion(2.0,3.0,4.0,5.0)
    q2 = Quaternion(6.0,7.0,8.0,9.0)
    qr = q1 * q2
    m1 = MultiVector3f(v2m((2.0,3.0,4.0,5.0))...)
    m2 = MultiVector3f(v2m((6.0,7.0,8.0,9.0))...)
    mr = m1 * m2
    @test mr == v2m(q2v((qr.s, qr.v1,qr.v2,qr.v3)))
end

# Octonions use a representation that is a particular permutation of
# our representation; I haven't figured out yet which it is.

# for j in 1:8, i in 1:8
#     const p = [1,2,3,4,5,7,8,6]
#     const f = [+1,+1,+1,+1,+1,+1,+1,+1]
#     function perm(tup)
#         a = [tup...]
#         permute!(a, p)
#         ((f.*a)...)
#     end
#     function iperm(tup)
#         a = f.*[tup...]
#         ipermute!(a, p)
#         (a...)
#     end
#     v1 = ([Float64(d==i) for d in 1:8]...)
#     v2 = ([Float64(d==j) for d in 1:8]...)
#     o1 = Octonion(perm(v1)...)
#     o2 = Octonion(perm(v2)...)
#     or = o1 * o2
#     ort = iperm([or.s, or.v1, or.v2, or.v3, or.v4, or.v5, or.v6, or.v7])
#     m1 = MultiVector3f(v1...)
#     m2 = MultiVector3f(v2...)
#     mr = m1 * m2
#     @show i j v1 v2 ort (mr...)
#     @test mr == MultiVector3f(ort...)
# end
# 
# let
#     o = (Octonion(2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0) *
#          Octonion(12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0))
#     @test (MultiVector3f(2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0) *
#            MultiVector3f(12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0) ==
#            MultiVector3f(o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7))
# end

@test (MultiVector3i(2, 0,0,0, 0, 0,0,0) *
       MultiVector3i(3, 0,0,0, 0, 0,0,0) ==
       MultiVector3i(6, 0,0,0, 0, 0,0,0))

@test (MultiVector3i(2, 0,0,0, 0, 0,0,0) *
       MultiVector3i(0, 3,4,0, 5, 0,0,0) ==
       MultiVector3i(0, 6,8,0, 10, 0,0,0))

@test (MultiVector3i(2, 0,0,0, 0, 0,0,0) *
       MultiVector3i(0, 0,0,3, 0, 4,5,0) ==
       MultiVector3i(0, 0,0,6, 0, 8,10,0))

@test (MultiVector3i(2, 0,0,0, 0,0,0, 0) *
       MultiVector3i(0, 0,0,0, 0,0,0, 3) ==
       MultiVector3i(0, 0,0,0, 0,0,0, 6))

@test (MultiVector3i(0, 2,3,0, 4,0,0, 0) *
       MultiVector3i(5, 0,0,0, 0,0,0, 0) ==
       MultiVector3i(0, 10,15,0, 20,0,0, 0))

@test (MultiVector3i(0, 2,3,0, 4,0,0, 0) *
       MultiVector3i(0, 5,6,0, 7,0,0, 0) ==
       MultiVector3i(-2*5-3*6-4*7, 0,0,2*6-3*5, 0,2*7-4*5,3*7-4*6, 0))

@test (MultiVector3i(0, 2,3,0, 4,0,0, 0) *
       MultiVector3i(0, 0,0,5, 0,6,7, 0) ==
       MultiVector3i(0, 3*5+4*6,-2*5+4*7,0, -2*6-3*7,0,0, 2*7-3*6+4*5))

@test (MultiVector3i(0, 2,3,0, 4,0,0, 0) *
       MultiVector3i(0, 0,0,0, 0,0,0, 5) ==
       MultiVector3i(0, 0,0,-4*5, 0,3*5,-2*5, 0))

@test (MultiVector3i(0, 0,0,2, 0,3,4, 0) *
       MultiVector3i(5, 0,0,0, 0,0,0, 0) ==
       MultiVector3i(0, 0,0,10, 0,15,20, 0))

@test (MultiVector3i(0, 0,0,2, 0,3,4, 0) *
       MultiVector3i(0, 5,6,0, 7,0,0, 0) ==
       MultiVector3i(0, -2*6-3*7,2*5-4*7,0, 3*5+4*6,0,0, 2*7-3*6+4*5))

@test (MultiVector3i(0, 0,0,2, 0,3,4, 0) *
       MultiVector3i(0, 0,0,5, 0,6,7, 0) ==
       MultiVector3i(-2*5-3*6-4*7, 0,0,3*7-4*6, 0,-2*7+4*5,2*6-3*5, 0))

@test (MultiVector3i(0, 0,0,2, 0,3,4, 0) *
       MultiVector3i(0, 0,0,0, 0,0,0, 5) ==
       MultiVector3i(0, -4*5,3*5,0, -2*5,0,0, 0))

@test (MultiVector3i(0, 0,0,0, 0,0,0, 2) *
       MultiVector3i(3, 0,0,0, 0,0,0, 0) ==
       MultiVector3i(0, 0,0,0, 0,0,0, 6))

@test (MultiVector3i(0, 0,0,0, 0,0,0, 2) *
       MultiVector3i(0, 3,4,0, 5,0,0, 0) ==
       MultiVector3i(0, 0,0,-2*5, 0,2*4,-2*3, 0))

@test (MultiVector3i(0, 0,0,0, 0,0,0, 2) *
       MultiVector3i(0, 0,0,3, 0,4,5, 0) ==
       MultiVector3i(0, -2*5,2*4,0, -2*3,0,0, 0))

@test (MultiVector3i(0, 0,0,0, 0,0,0, 2) *
       MultiVector3i(0, 0,0,0, 0,0,0, 3) ==
       MultiVector3i(6, 0,0,0, 0,0,0, 0))
