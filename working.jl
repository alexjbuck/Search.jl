using Pkg
Pkg.activate("Search")
using Search
using Distributions
using LinearAlgebra
using Plots
initDist = MvNormal(Array{Float64}([1000, 1000, 100, 100]))

searchDist = MvNormal((Array{Float64}([1000,1000])))

searchSchedule = [0,1];
function update!(targets,Δt)
    targets[:,1:2] = targets[:,1:2] .+ targets[:,3:4] .* Δt
end
Search.initializeScenario(initDist,searchDist,searchSchedule,Search.linearTargetModel!)
N = 10;
targets = Array{Float64}(undef,(N,4))
for n = 1:N
    targets[n,:] = rand(initDist)
end
lim = 10000
xrange = (-lim,lim);
yrange = (-lim,lim);
scatter(targets[:,1],targets[:,2];xlim=xrange,ylim=yrange)


function initTargets(initDist::MvNormal;N=100)
    targets = Array{Float64}(undef,(N,4))
    for n = 1:N
        targets[n,:] = rand(initDist)
    end
    targets = [targets zeros(N,1)]
    return targets
end

function isFound(x1,x2,threshold)
    norm(x2.-x1)<=threshold ? true : false
end

Δt = 1/1;
targets = initTargets(initDist;N);
scatter(targets[:,1],targets[:,2];xlim=xrange,ylim=yrange)
for t = 1:Δt:50
    update!(targets,Δt)
    p = scatter(targets[:,1],targets[:,2];xlim=xrange,ylim=yrange)
    for i in 1:N
        if isFound(targets[i,1:2],[0 0],10000)
            targets[i,:] .= NaN;
            targets[i,end] = 1;
        end
    end
    if all(i -> i==1,targets[:,5])
        println("All targets found")
        break
    end
    p = scatter(targets[:,1],targets[:,2];xlim=xrange,ylim=yrange)
    display(p)
    sleep(.1)
end


dims(x) = length(size(x))

function gridDist(d::MvNormal;δ = .1)
    Δ = 2*sqrt(det(cov(d)))
    i = -Δ:δ:Δ;
    j = -Δ:δ:Δ;
    z = [pdf(d,[i,j]) for i in i, j in j]
    return i,j,z
end

function addGrid(i,j,A::AbstractArray,k,l,B::AbstractArray,Δ)
    x = minimum([i...,k...]):Δ :maximum([i...,k...]);
    y = minimum([j...,l...]):Δ :maximum([j...,l...]);
    _,_,pad_A = padGrid(x,y,i,j,A,Δ)
    _,_,pad_B = padGrid(x,y,k,l,B,Δ)
    G = pad_A .+ pad_B
    return x,y,G
end

function multGrid(i,j,A::AbstractArray,k,l,B::AbstractArray,Δ)
    x,y,overlap_A = overlapGrid(k,l,i,j,A,Δ)
    _,_,overlap_B = overlapGrid(i,j,k,l,B,Δ)
    G = overlap_A .* overlap_B
    return x,y,G
end

function padGrid(x,y,i,j,A,Δ)
    a = abs(minimum(i) - minimum(x))/Δ 
    b = length(y)
    pad_upper = zeros(convert.(Int,[a,b])...)
    println(typeof(a))
    println(typeof(b))
    a = abs(maximum(i) - maximum(x))/Δ 
    pad_lower = zeros(convert.(Int,[a,b])...)
    println(typeof(a))
    println(typeof(b))
    a = length(i)
    b = abs(minimum(j) - minimum(y))/Δ 
    pad_left = zeros(convert.(Int,[a,b])...)
    println(typeof(a))
    println(typeof(b))
    b = abs(maximum(j) - maximum(y))/Δ 
    pad_right = zeros(convert.(Int,[a,b])...)
    println(typeof(a))
    println(typeof(b))
    pad_A = [pad_upper;pad_left A pad_right;pad_lower]
    return x,y,pad_A
end

function overlapGrid(k,l,i,j,A,Δ)
    x = maximum( [minimum(k),minimum(i)] ):Δ:minimum( [maximum(k),maximum(i)] )
    y = maximum( [minimum(l),minimum(j)] ):Δ:minimum( [maximum(l),maximum(j)] )
    xi = indexin([minimum(x),maximum(x)],i)
    yi = indexin([minimum(y),maximum(y)],j)
    overlap_A = A[xi[1]:xi[2],yi[1]:yi[2]]
    return x,y,overlap_A
end