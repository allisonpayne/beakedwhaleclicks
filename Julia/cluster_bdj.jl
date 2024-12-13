using Feather
using DataFrames
using FileIO
using Statistics, StatsBase
using Plots
using StatsPlots
using DSP
using Clustering
using Random




dt = 1/256000

features_root = "C:\\Users\\bjones\\Desktop\\MarsClickProcessing\\features"
println("Listing detected feature files...")
featurefiles_vec = []
for (root, dirs, files) in walkdir(features_root)
    ff = filter(f -> endswith(f, ".feather"), files)
    append!(featurefiles_vec, joinpath.(features_root, root, ff))
end

featureslist = []

for f in featurefiles_vec
    println(f)
    featurestemp=Feather.read(f)
    push!(featureslist,featurestemp)
end

features = vcat(featureslist...)

features.logduration = log.(max.(dt, features.duration .+ 2e-6randn(size(features, 1))))

feature_names = [:logduration, :peak_freq, :peak2_freq, :notch1_freq, :notch2_freq,
    :peak_bw, :freq05, :freq50, :freq95, :mean_freq]
X = Matrix(features[:, feature_names]) # Clustering functions expect a Matrix, not a DataFrame
nfeat = size(X, 2)

tiles1 = []
for i in 1:nfeat
    x = X[:, i]
    limits = quantile(real.(x), [0.0001, 0.9999])
    p = density(x, yticks=false, legend=false,
         title=feature_names[i])
    push!(tiles1, p)
end

p = plot(tiles1..., layout=(2, 5), size=(1400, 1000))

feat_i = 1:nfeat
tiles2 = []
for i in feat_i, j in feat_i
    if i == j
        p = density(real.(X[:, i]), yticks=false, legend=false, title=feature_names[i])
    else
        p = histogram2d(real.(X[:, j]), real.(X[:, i]), nbins=100,
            yticks=false, xticks=false, legend=false)
    end
    push!(tiles2, p)
end
p = plot(tiles2..., layout=(nfeat, nfeat), size=(3000, 3000))
savefig(p, "graphics\\feature_histmatrix.png")

###############################################################################
Random.seed!(1052544878090191442)
min_x = minimum(X, dims=1)
Xnorm = X .- min_x
max_x = maximum(Xnorm, dims=1)
Xnorm = Xnorm ./  max_x
#feat_i = [1,2,3,4,5,6,7,8,9,10] # Just do clustering with four of the features
feat_i = [1,2,6,8] # Just do clustering with four of the features
nsub = 22_200 # Reduce the size of the sample to make things run a little faster
sub_i = rand(1:size(Xnorm,1), nsub)
Xnorm_sub = Xnorm[sub_i, feat_i]
Xnorm_sub = collect(Xnorm_sub')
Waveform_sub = waveforms[sub_i,:]

nclust = 15
clust = fuzzy_cmeans(Xnorm_sub, nclust, 2, display=:iter, maxiter=200)

realcenters = clust.centers .* max_x[feat_i] .+ min_x[feat_i]
tiles3 = []
n = length(feat_i)
for i in 1:n, j in 1:n
    if i == j
        p = density(X[:, feat_i[i]], yticks=false, legend=false, title=feature_names[feat_i[i]])
    else
        p = histogram2d(X[:, feat_i[j]], X[:, feat_i[i]], nbins=100,
            yticks=false, legend=false)
        scatter!(p, realcenters[j,:], realcenters[i,:],
            series_annotations=text.(1:nclust, color=:white))
    end
    push!(tiles3, p)
end
p = plot(tiles3..., layout=(length(feat_i), length(feat_i)), size=(3000, 3000));


