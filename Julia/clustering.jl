using JLD2, FileIO
using CSV
using DataFrames
using Statistics
using DSP
using StableRNGs
using UMAP
using Clustering
using ConcaveHull
using Alert
using StatsPlots

include(joinpath(@__DIR__, "deimos_detector.jl"))

data_dir = "/media/sam/Sam01/"

# Waveforms
fs = 256e3
dt = 1/fs
waveforms = collect(load(joinpath(data_dir, "click_subsample.jld2"), "waveforms")')
waveforms .-= mean(waveforms, dims=1)
s = vec(std(waveforms, dims=1))
histogram(log10.(s))
waveforms ./= s'
# waveforms = waveforms[:, s .> quantile(s, 0.25)]

# Calculating spectrum of each click
n = 256
x = waveforms[:, 1]
x = [zeros(256); x; zeros(256)]
pg = welch_pgram(x, n, fs=fs)
pg.freq
plot(pg.freq/1e3, pg.power)

spectra = zeros(length(pg.power), size(waveforms, 2))
for i in 1:size(waveforms, 2)
    x = waveforms[:, i]
    x = [zeros(256); x; zeros(256)]
    pg = welch_pgram(x, n, fs=fs)
    spectra[:, i] .= pg.power
end

StableRNGs.seed!(123)
idx_deimos = findall(s -> is_deimos(pg.freq, s), eachcol(spectra))
idx_1 = setdiff(1:size(waveforms, 2), idx_deimos)
idx_1 = sample(idx_1, 50_000, replace=false)
waveforms1 = waveforms[:, idx_1]
spectra1 = spectra[:, idx_1]

rand_spectra_plots = map(rand(axes(spectra1, 2), 24)) do i
    plot(pg.freq/1e3, spectra1[:, i], label="")
end
plot(rand_spectra_plots..., xlabel="", ylabel="", yticks=false, 
    size=(1200, 700), layout=(6, 4), dpi=150)
savefig(joinpath(@__DIR__, "../graphics/rand_spectra.png"))


#=
Distance/feature engineering section
=# 

include(joinpath(@__DIR__, "distance_metrics.jl"))
include(joinpath(@__DIR__, "rissos_pwsd.jl"))

i_rissos = findall(s -> rissos(pg.freq, s), eachcol(spectra1))
i_pwsd = findall(s -> pwsd(pg.freq, s), eachcol(spectra1))
i_other = setdiff(1:size(spectra1, 2), (union(i_rissos, i_pwsd)))

s_rissos = vec(mean(spectra1[:, i_rissos], dims=2))
s_pwsd = vec(mean(spectra1[:, i_pwsd], dims=2))

plot([s_rissos, s_pwsd], labels=["Risso's" "PWSD"], color=[2 3], xlim=(0, 80),
    xlabel="Frequency (kHz)", yticks=false, linewidth=3, size=(400, 300), dpi=300)
savefig(joinpath(@__DIR__, "../graphics/rissos_v_pwsd.png"))

g = zeros(size(spectra1, 2))
g[i_rissos] .= 2
g[i_pwsd] .= 3

distances = [
    Euclidean(),
    SqEuclidean(),
    ExpEuclidean(0.5),
    Minkowski(0.25),
    CosineDist(),
    BrayCurtis(),
    EarthMovers(),
    EarthMovers(2),
    CorrDist(),
    CrossCorrDist()
]



dist_matrix = zeros(3, 3, length(distances))
nreps = 10_000
for i in 1:nreps
    for (j, dist) in enumerate(distances)
        sr1 = spectra1[:, rand(i_rissos)]
        sr2 = spectra1[:, rand(i_rissos)]
        sp1 = spectra1[:, rand(i_pwsd)]
        sp2 = spectra1[:, rand(i_pwsd)]
        so1 = spectra1[:, rand(i_other)]
        so2 = spectra1[:, rand(i_other)]
        dist_matrix[1, 1, j] += dist(sr1, sr2)
        dist_matrix[2, 2, j] += dist(sp1, sp2)
        dist_matrix[3, 3, j] += dist(so1, so2)
        dist_matrix[2, 1, j] += dist(sp1, sr1)
        dist_matrix[3, 1, j] += dist(so1, sr1)
        dist_matrix[3, 2, j] += dist(so1, sp1)
    end
end
dist_matrix ./= nreps
p_dists = map(enumerate(distances)) do (i, dist)
    heatmap(dist_matrix[:, :, i], yflip=true, title=string(dist))
end
plot(p_dists..., size=(800, 600))

otherdists = map(distances) do d
    Δ = map(1:10_000) do _
        s1 = spectra1[: , rand(i_other)]
        s2 = spectra1[:, rand(i_other)]
        Distances.evaluate(d, s1, s2)
    end
    return Δ
end

dolphindists = map(distances) do d
    Δ = map(1:10_000) do _
        s1 = spectra1[:, rand(i_rissos)]
        s2 = spectra1[:, rand(i_pwsd)]
        Distances.evaluate(d, s1, s2)
    end
    return Δ
end

distplots = map(enumerate(distances)) do (i, dist)
    p = density(otherdists[i], title=string(dist), legend=false)
    density!(p, dolphindists[i])
    p
end
plot(distplots..., size=(800, 800))

reldists = mean.(dolphindists) ./ mean.(otherdists)

spec_weights = vec(mean(spectra1, dims=2))
spec_weights = spec_weights ./ sum(spec_weights)
plot(spec_weights)
# dolphin_freqs = [19:23; 25:28; 30:31; 33] * 1e3
# dolphin_freqs = [21:23; 28:31] * 1e3
# spec_weights = ones(length(pg.freq))
# spec_weights[in(dolphin_freqs).(pg.freq)] .*= 5
# spec_weights[in((30:31)*1e3).(pg.freq)] .*= 5
# spec_weights = spec_weights ./ sum(spec_weights)

spec_weights1 = abs.(s_rissos .- s_pwsd)
# spec_weights .+= 1/129
spec_weights1 ./= sum(spec_weights1)
plot(spec_weights1)
spec_weights = spec_weights .+ spec_weights1
spec_weights = spec_weights ./ sum(spec_weights)
plot(spec_weights)


StableRNGs.seed!(123)
embedding = umap(waveforms1, 2, n_neighbors=6,
    # metric=WeightedExpEuclidean(spec_weights, 0.5),
    metric=CorrDist(),
    min_dist=0.1, init=:random)
alert("UMAP done.")
scatter(embedding[1, :], embedding[2, :], group=g, label=["Other" "Rissos" "PWSD"],
    markersize=1, markerstrokewidth=0, size=(800, 800))

# Parameter search for UMAP
n_neighbors = 3:8
# min_dist = [0.05, 0.1, 0.125, 0.15, 0.2, 0.25, 0.5]
# min_dist = range(0.001, 0.2, length=6)
min_dist = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
params = [(n_neighbors=n, min_dist=d) for n in n_neighbors, d in min_dist]

embeddings = map(params) do p
    println(p)
    StableRNGs.seed!(123)
    umap(waveforms1, 2; p..., 
        metric=CorrDist(),
        init=:random)
end

pp = map(enumerate(params)) do (i, p)
    embed = embeddings[i]
    title = "n = $(p.n_neighbors), d = $(p.min_dist)"
    p = scatter(embed[1, :], embed[2, :], size=(800, 800), xticks=false, yticks=false,
        markersize=0.5, markerstrokewidth=0, label="", title=title)
    scatter!(p, embed[1, i_pwsd], embed[2, i_pwsd], label="",
        markerstrokewidth=0, markersize=1, color=3) 
    scatter!(p, embed[1, i_rissos], embed[2, i_rissos], label="",
        markerstrokewidth=0, markersize=1, color=2) 
end
p = plot(pp..., layout=(6, 6), size=(1400, 1400));
savefig(p, joinpath(@__DIR__, "../graphics/umap_parameters_q25+_waveforms_CorrDist().png"))

######################################################

StableRNGs.seed!(123)
umap_model = UMAP_(spectra1, 2, n_neighbors=6, 
    metric=EarthMovers(),
    min_dist=0.001, init=:random)
alert("UMAP done.")
save(joinpath(@__DIR__, "../data/umap_model.jld2"), Dict("umap_model" => umap_model))
umap_model = load(joinpath(@__DIR__, "../data/umap_model.jld2"), "umap_model")
embedding = umap_model.embedding


scatter(embedding[1, :], embedding[2, :], markersize=1, markerstrokewidth=0,
    legend=false, size=(800, 800), dpi=150)
savefig(joinpath(@__DIR__, "../graphics/umap_embedding.png"))

scatter(embedding[1, :], embedding[2, :], group=g, label=["Other" "Rissos" "PWSD"],
    markersize=1, markerstrokewidth=0, size=(800, 800))
savefig(joinpath(@__DIR__, "../graphics/umap_embedding_rissos_pwsd.png"))
# savefig(joinpath(@__DIR__, "../graphics/umap_embedding_rissos_pwsd.png"))


#= 
Density-based clustering of points
=#
clusters = dbscan(embedding, 0.5, min_neighbors=15)#, min_cluster_size=10)
nclusters(clusters)

all_counts = [sum(clusters.assignments .== 0); clusters.counts]
cluster_labels = DataFrame(
    assignments = 0:nclusters(clusters),
    label = invperm(sortperm(all_counts, rev=true))
)

results = DataFrame(embedding' , [:e1, :e2])
results.assignments = clusters.assignments

results = leftjoin(results, cluster_labels, on=:assignments)
CSV.write(joinpath(@__DIR__, "../data/clusters.csv"), results)

centroids = combine(
    groupby(results, [:assignments, :label]),
    :e1 => mean, :e2 => mean
)
centroids = subset(centroids, :assignments => x -> x .> 0)

cluster_hulls = map(1:nclusters(clusters)) do i
    points = eachcol(embedding)[i .== clusters.assignments]
    concave_hull(points, 100)
end
cluster_hulls = cluster_hulls[sortperm(clusters.counts, rev=true)]
save(joinpath(@__DIR__, "../data/cluster_hulls.jld2"), Dict("cluster_hulls" => cluster_hulls))

p_clust = scatter(embedding[1, :], embedding[2, :], group=results.label,
    markersize=0.5, markerstrokewidth=0, legend=false, size=(800, 800));
for (i, h) in enumerate(cluster_hulls)
    plot!(p_clust, h, color=i)
end
scatter!(p_clust, centroids.e1_mean .+ 0.75, centroids.e2_mean, series_annotations=centroids.label,
    markersize=0, markerstrokewidth=0, color=:white);
p_clust
savefig(p_clust, joinpath(@__DIR__, "../graphics/clusters1.png"))

#= 
Average spectra and waveforms
=#
spectra2 = spectra[:, setdiff(1:size(spectra, 2), idx_deimos)]

avg_spectra = map(1:nclusters(clusters)) do i
    s = mean(spectra1[:, findall(i .== results.label)], dims=2)
    s[47 .< pg.freq/1e3 .< 53] .= NaN
    s
end
plot(pg.freq/1e3, avg_spectra, label=[1:nclusters(clusters);]', color=[1:nclusters(clusters);]',
    xticks=0:10:120,
    layout=(3, 5), yticks=false, xlims=(0, 100), size=(1200, 700), dpi=150)
savefig(joinpath(@__DIR__, "../graphics/average_spectra.png"))


function maxcorr_mean(X, lags=-100:100)
    n = size(X, 2)
    shifts = zeros(Int, n)
    signs = ones(n)
    sum = zero(X[:, 1])
    for i in 2:n
        cc = crosscor(X[:, 1], X[:, i])
        abscc = abs.(cc)
        imax = argmax(abscc)
        shifts[i] = lags[imax]
        signs[i] = sign(cc[imax])
        sum .+= circshift(X[:, i], -shifts[i]) .* signs[i]
    end
    return sum / n
end

# using MultivariateStats
# i = 9
# W = waveforms1[:, findall(i .== results.label)]
# plot(W)
# pca = fit(PCA, W, pratio=0.95);
# Φ = projection(pca)
# A = predict(pca, W)
# W_smooth = reconstruct(pca, A)
# # W_smooth = Φ * mean(A, dims=2)
# plot(maxcorr_mean(W))
# plot!(mean(W_smooth, dims=2))

avg_waveforms = map(1:nclusters(clusters)) do i
    maxcorr_mean(waveforms1[:, findall(i .== results.label)])
    # W = waveforms1[:, findall(i .== results.label)]
end
tt = (1:512) / 256e3 * 1e3
plot(tt, avg_waveforms, yticks=false, label=[1:nclusters(clusters);]', color=[1:nclusters(clusters);]',
    layout=(3, 5), size=(1200, 700))
savefig(joinpath(@__DIR__, "../graphics/average_waveforms.png"))

spectra_table = map(enumerate(avg_spectra)) do (i, spec)
    df = DataFrame(cluster = i, freq = pg.freq, power = vec(spec))
    fmax = df.freq[argmax(replace(spec, NaN => -Inf))] / 1e3
    p = plot(df.freq / 1e3, df.power, label="Peak = $fmax kHz", xticks=0:10:120,
        size=(1000, 600))
    vline!(p, [22.13, 26.75, 33.2, 37.29], linestyle=:auto, label="PWSD")
    vline!(p, [22.05, 25.58, 30.34, 39.04], linestyle=:auto, label="Risso's")
    vline!(p, [16.4], linestyle=:auto, label="Baird's")
    vline!(p, [40.2], linestyle=:auto, label="Cuvier's")
    vline!(p, [34.4], linestyle=:auto, label="Blainville's")
    ipad = lpad(string(i), 2, "0")
    savefig(p, joinpath(@__DIR__, "../graphics/spectra/spectrum_$ipad.png"))
    return df
end;
spectra_table = vcat(spectra_table...)
CSV.write(joinpath(@__DIR__, "../data/average_spectra.csv"), spectra_table)

waveform_table = map(enumerate(avg_waveforms)) do (i, wav)
    DataFrame(cluster = i, t = (1:length(wav)) ./ fs, waveform = vec(wav))
end
waveform_table = vcat(waveform_table...)
CSV.write(joinpath(@__DIR__, "../data/average_waveforms.csv"), spectra_table)
