module GeometricAlgebra

using FlexibleArrays

# function rank_lengths(dim::Int)
#     lengths = Int[]
#     len = 1
#     for i in 0:dim
#         push!(lengths, len)
#         @assert len * (dim-i) % (i+1) == 0
#         len = len * (dim-i) ÷ (i+1)
#         @assert len >= 0
#         @assert (len == 0) == (i == dim)
#     end
#     @assert length(lengths) == dim+1
#     @assert sum(lengths) == 2^dim
#     lengths
# end



export AbstractMultiVector
abstract AbstractMultiVector{T,D}

export MultiVector
immutable MultiVectorImpl{D,M,F,T,CT} <: AbstractMultiVector{T,D}
    elems::CT
    function MultiVectorImpl()
        bnds = [(:) for d in 1:D]
        new(CT(bnds...))
    end
    function MultiVectorImpl(elems::CT)
        new(elems)
    end
    function MultiVectorImpl(elems::T...)
        new(CT(nothing, elems))
    end
    # function MultiVectorImpl(elems::Tuple)
    #     @assert length(elems) == D+1
    #     for d in 0:D
    #         @assert (length[elems[d+1]] ==
    #                  factorial(D) ÷ (factorial(D-d) * factorial(d)))
    #         @assert eltype(elems[d+1]) === Int
    #     end
    #     counts = [0 for d in 0:D]
    #     ...
    #     new(CT(nothing, elems))
    # end
end

function MultiVector{D,DD,T}(d::Int,
                             Moduli::NTuple{D,Int},
                             Factors::NTuple{DD,Int},
                             ::Type{T})
    @assert d == D
    @assert DD == prod(Moduli)^3
    C = ImmutableArray([0:m-1 for m in Moduli]...)
    MultiVectorImpl{D, Moduli, Factors, T, C{T}}
end

function MultiVector{D,T}(d::Int,
                          Moduli::NTuple{D,Int},
                          Parities::NTuple{D,Int},
                          Commutators::NTuple{D,NTuple{D,Int}},
                          ::Type{T})
    @assert d == D
    for e in 1:D, d in 1:D
        @assert Commutators[d][e] == Commutators[e][d]
    end
    Factors =
        FlexArray([0:m-1 for m in Moduli]...,
                  [0:m-1 for m in Moduli]...,
                  [0:m-1 for m in Moduli]...){Int}([(:) for m in Moduli]...,
                                                   [(:) for m in Moduli]...,
                                                   [(:) for m in Moduli]...)
    for ind in eachindex(Factors)
        Factors[ind] = 0
    end
    lind = CartesianIndex([0 for m in Moduli]...)
    uind = CartesianIndex([m-1 for m in Moduli]...)
    for i in CartesianRange(lind, uind)
        for j in CartesianRange(lind, uind)
            k = i
            f = +1
            for jd in 1:D
                for kd in D:-1:jd+1
                    f *= (Commutators[kd][jd] ^ k[kd]) ^ j[jd]
                end
                kval = k[jd] + j[jd]
                if kval >= Moduli[jd]
                    kval -= Moduli[jd]
                    f *= Parities[jd]
                end
                k = CartesianIndex(setindex(k.I, kval, jd))
            end

            symbols = Int[]
            for d in 1:D
                for n in 1:i[d]
                    push!(symbols, d)
                end
            end
            for d in 1:D
                for n in 1:j[d]
                    push!(symbols, d)
                end
            end
            f2 = +1
            while true
                issorted = true
                for n in 2:length(symbols)
                    if symbols[n] < symbols[n-1]
                        s1 = symbols[n-1]
                        s2 = symbols[n]
                        symbols[n-1] = s2
                        symbols[n] = s1
                        f2 *= Commutators[s1][s2]
                        issorted = false
                        break
                    end
                end
                issorted && break
            end
            while true
                nodups = true
                for n in 2:length(symbols)
                    if symbols[n] == symbols[n-1]
                        s1 = symbols[n]
                        symbols = vcat(symbols[1:n-2], symbols[n+1:end])
                        f2 *= Parities[s1]
                        nodups = false
                        break
                    end
                end
                nodups && break
            end
            for d in 1:D
                @assert count(x->x==d, symbols) == k[d]
            end
            for d in 1:D
                @assert mod(i[d] + j[d] - k[d], Moduli[d]) == 0
            end            
            @assert f2 == f

            Factors[i.I..., j.I..., k.I...] += f
        end
    end
    MultiVector(D, Moduli, (Factors...), T)
end

function MultiVector{T}(D::Int, ::Type{T})
    Moduli = FlexArray(1:D){Int}(:)
    for d in 1:D
        Moduli[d] = 2
    end
    Parities = FlexArray(1:D){Int}(:)
    for d in 1:D
        Parities[d] = -1
    end
    Commutators = FlexArray(1:D,1:D){Int}(:,:)
    for e in 1:D, d in 1:D
        Commutators[d,e] = d==e ? 0 : -1
    end
    MultiVector(D,
                (Moduli...),
                (Parities...),
                ntuple(e -> ntuple(d -> Commutators[d,e], D), D),
                T)
end



import Base: ==, +, -, *
import Base: done, eachindex, eltype, fill, getindex, length, ndims, next, similar, size, start, zero
import FlexibleArrays: setindex
export setindex

eltype{D,M,F,T,CT}(::Type{MultiVectorImpl{D,M,F,T,CT}}) = eltype(CT)
ndims{D,M,F,T,CT}(::Type{MultiVectorImpl{D,M,F,T,CT}}) = ndims(CT)
size{D,M,F,T,CT}(::Type{MultiVectorImpl{D,M,F,T,CT}}) = size(CT)
length{D,M,F,T,CT}(::Type{MultiVectorImpl{D,M,F,T,CT}}) = length(CT)
 
ndims(m::AbstractMultiVector) = ndims(typeof(m))
size(m::AbstractMultiVector) = size(typeof(m))
eltype(m::AbstractMultiVector) = eltype(typeof(m))
length(m::AbstractMultiVector) = length(typeof(m))

eachindex(m::AbstractMultiVector) = eachindex(m.elems)

start(m::AbstractMultiVector) = start(m.elems)
done(m::AbstractMultiVector, state) = done(m.elems, state)
next(m::AbstractMultiVector, state) = next(m.elems, state)
 
 
 
zero{MV <: AbstractMultiVector}(::Type{MV}) = MV()
similar(m::AbstractMultiVector) = zero(typeof(m))



function getindex(m::AbstractMultiVector, inds::Int...)
    @assert length(inds) == ndims(m)
    m.elems[inds...]
end
function getindex{D}(m::AbstractMultiVector, inds::CartesianIndex{D})
    m[ntuple(d->inds[d], D)...]
end

function setindex(m::AbstractMultiVector, val, inds::Int...)
    @assert length(inds) == ndims(m)
    typeof(m)(setindex(m.elems, val, inds...))
end
function setindex{D}(m::AbstractMultiVector, val, inds::CartesianIndex{D})
    setindex(m, val, ntuple(d->inds[d], D)...)
end



function (==){T,D}(m1::AbstractMultiVector{T,D}, m2::AbstractMultiVector{T,D})
    r = true
    for i in eachindex(m1)
        r &= m1[i] == m2[i]
    end
    r
end

function (+)(m1::AbstractMultiVector)
    r = similar(m1)
    for i in eachindex(r)
        r = setindex(r, + m1[i], i)
    end
    r
end

function (+){T,D}(m1::AbstractMultiVector{T,D}, m2::AbstractMultiVector{T,D})
    r = similar(m1)
    for i in eachindex(r)
        r = setindex(r, m1[i] + m2[i], i)
    end
    r
end

function (-)(m1::AbstractMultiVector)
    r = similar(m1)
    for i in eachindex(r)
        r = setindex(r, - m1[i], i)
    end
    r
end

function (-){T,D}(m1::AbstractMultiVector{T,D}, m2::AbstractMultiVector{T,D})
    r = similar(m1)
    for i in eachindex(r)
        r = setindex(r, m1[i] - m2[i], i)
    end
    r
end

# TODO: Use promotion instead of accepting Number
function (*)(x::Number, m1::AbstractMultiVector)
    r = similar(m1)
    for i in eachindex(r)
        r = setindex(r, x * m1[i], i)
    end
    r
end

function (*)(m1::AbstractMultiVector, x::Number)
    r = similar(m1)
    for i in eachindex(r)
        r = setindex(r, m1[i] * x, i)
    end
    r
end

function (*){D,M,F,T,CT}(m1::MultiVectorImpl{D,M,F,T,CT},
                         m2::MultiVectorImpl{D,M,F,T,CT})
    f = ImmutableArray([0:m-1 for m in M]...,
                       [0:m-1 for m in M]...,
                       [0:m-1 for m in M]...){Int}(nothing, F)
    r = MultiVectorImpl{D,M,F,T,CT}()
    lind = CartesianIndex([0 for m in M]...)
    uind = CartesianIndex([m-1 for m in M]...)
    for k in CartesianRange(lind, uind)
        s = zero(T)
        for j in CartesianRange(lind, uind)
            for i in CartesianRange(lind, uind)
                if f[i,j,k] != 0
                    s += f[i,j,k] * m1[i] * m2[j]
                end
            end
        end
        r = setindex(r, s, k)
    end
    r
end

end
