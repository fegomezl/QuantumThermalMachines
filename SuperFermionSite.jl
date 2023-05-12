using ITensors

#Create the SuperFermion site with the corresponding local dimension
ITensors.space(::SiteType"SuperFermion") = 4

#Create alias to the possible states
ITensors.state(::StateName"00", ::SiteType"SuperFermion") = [1., 0., 0., 0.]
ITensors.state(::StateName"01", ::SiteType"SuperFermion") = [0., 1., 0., 0.]
ITensors.state(::StateName"10", ::SiteType"SuperFermion") = [0., 0., 1., 0.]
ITensors.state(::StateName"11", ::SiteType"SuperFermion") = [0., 0., 0., 1.]

ITensors.state(::StateName"Vacuum", ::SiteType"SuperFermion") = [1., 0., 0., 1.]
ITensors.state(::StateName"NormalizedVacuum", ::SiteType"SuperFermion") = [0.5, 0., 0., 0.5]

#Create operators
#P:Physical A:Ancilla

#Local operators
function ITensors.op(::OpName"I", ::SiteType"SuperFermion", s::Index)
    mat = [1. 0. 0. 0.
           0. 1. 0. 0.
           0. 0. 1. 0.
           0. 0. 0. 1.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"bᴾ", ::SiteType"SuperFermion", s::Index)
    mat = [0. 0. 1. 0.
           0. 0. 0. 1.
           0. 0. 0. 0.
           0. 0. 0. 0.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"bᴬ", ::SiteType"SuperFermion", s::Index)
    mat = [0. 1. 0. 0.
           0. 0. 0. 0.
           0. 0. 0. 1.
           0. 0. 0. 0.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"b†ᴾ", ::SiteType"SuperFermion", s::Index)
    mat = [0. 0. 0. 0.
           0. 0. 0. 0.
           1. 0. 0. 0.
           0. 1. 0. 0.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"b†ᴬ", ::SiteType"SuperFermion", s::Index)
    mat = [0. 0. 0. 0.
           1. 0. 0. 0.
           0. 0. 0. 0.
           0. 0. 1. 0.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"nᴾ", ::SiteType"SuperFermion", s::Index)
    mat = [0. 0. 0. 0.
           0. 0. 0. 0.
           0. 0. 1. 0.
           0. 0. 0. 1.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"nᴬ", ::SiteType"SuperFermion", s::Index)
    mat = [0. 0. 0. 0.
           0. 1. 0. 0.
           0. 0. 0. 0.
           0. 0. 0. 1.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"JWᴾ", ::SiteType"SuperFermion", s::Index)
    mat = [1. 0.  0.  0.
           0. 1.  0.  0.
           0. 0. -1.  0.
           0. 0.  0. -1.]
    return itensor(mat, s', s)
end

function ITensors.op(::OpName"JWᴬ", ::SiteType"SuperFermion", s::Index)
    mat = [1.  0. 0.  0.
           0. -1. 0.  0.
           0.  0. 1.  0.
           0.  0. 0. -1.]
    return itensor(mat, s', s)
end

# 2-site operators
function ITensors.op(::OpName"SWAP", ::SiteType"SuperFermion", s1::Index, s2::Index)
    mat = [1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.
           0.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.  0.
           0.  0.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.
           0.  0.  0.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.  0.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.  0.
           0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.]
    return itensor(mat, s2', s1', s2, s1)
end
