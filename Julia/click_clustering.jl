using FileIO
using Statistics, StatsBase
using DataFrames
using Plots, StatsPlots
using DSP
using Clustering
using Random
using FeatherFiles

dt = 1/256000

#waveforms = load("F:\\click_subsample.jld2", "waveforms")
#waveforms = waveforms ./ std(waveforms, dims=2)

features = Feather.Read("data/features_20190529.feather")
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
feat_i = [1,2,6,8] # Just do clustering with four of the features
nsub = 50_000 # Reduce the size of the sample to make things run a little faster
sub_i = rand(1:size(Xnorm,1), nsub)
Xnorm_sub = Xnorm[sub_i, feat_i]
Xnorm_sub = collect(Xnorm_sub')

# clusts1 = [fuzzy_cmeans(Xnorm_sub, nclust, 2, display=:iter, maxiter=200) for nclust in 2:20]
# clusts2 = [fuzzy_cmeans(Xnorm_sub, nclust, 2, display=:iter, maxiter=200) for nclust in 2:20]
# function r2(clust::FuzzyCMeansResult, X::AbstractArray)
#     ss = 0
#     for i in 1:size(X, 2)
#         assign = argmax(clust.weights[i, :])
#         ss += sum(abs2, X[:, i] .- clust.centers[:, assign])
#     end
#     μ = mean(X, dims=2)
#     ss0 = sum(abs2, X .- μ)
#     return 1 - ss/ss0
# end
#
# plot(2:20, [max(r2(clusts1[i], Xnorm_sub), r2(clusts2[i], Xnorm_sub)) for i in 1:19])

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
savefig(p, "graphics\\clusters.png")


tiles3 = []
n = length(feat_i)
for i in 1:n, j in 1:n
    if i == j
        p = density(real.(X[:, feat_i[i]]), yticks=false, legend=false, title=feature_names[i])
    else
        h = fit(Histogram, (X[:, feat_i[i]], X[:, feat_i[j]]), nbins=100)
        logcount = log10.(h.weights)
        logcount[logcount .< 0] .= NaN
        p = heatmap(midpoints(h.edges[2]), midpoints(h.edges[1]), logcount,
            legend=false)
        scatter!(p, realcenters[j,:], realcenters[i,:],
            series_annotations=text.(1:nclust, color=:black))
    end
    push!(tiles3, p)
end
p = plot(tiles3..., layout=(length(feat_i), length(feat_i)), size=(2000, 2000));
savefig(p, "graphics\\clusters_logcount.png")

# clust.weights is an nsub x nclust Matrix. Each row contains the probabilities
# that one click belongs to each cluster.
histogram(maximum(clust.weights, dims=2))

# Assign each click to a cluster if it's highest membership probability is greater
# than 2x the even-odds probability (which is 1/nclust).  If it doesn't fit in
# any of the clusters, give it label 0 for unassigned.
labels = [maximum(row) > 2/nclust ? argmax(row) : 0 for row in eachrow(clust.weights)]
histogram(labels)


histogram2d(features.peak_freq/1e3, features.logduration, xticks=0:2:100,
    size=(1400, 1000))
histogram2d(features.peak_freq/1e3, features.freq05, xticks=0:2:100,
    size=(1400, 1000))
histogram2d(features.peak_freq/1e3, features.freq05, xticks=0:2:100,
    size=(1400, 1000))

nlabels = length(unique(labels))

p = plot(size=(1500, 1500));
for i in 0:nlabels
    x = X[sub_i, :][labels .== i, :]
    scatter!(p, x[:, feat_i[1]], x[:, feat_i[2]], label=i, markersize=0.1,
        markeralpha=1, markerstrokealpha=0)
end
p
savefig(p, "graphics\\cluster_assignments.png")

mean_waveforms = [vec(mean(waveforms[sub_i, :][labels .== i, :], dims=1))
    for i in 0:nlabels-1]
p = plot(mean_waveforms, layout=(4,4), legend=false, size=(1600, 1200))
savefig(p, "graphics\\mean_waveforms.png")

mean_pgs = [periodogram(w, fs=256e3) for w in mean_waveforms]
p = plot(mean_pgs[1].freq/1e3, [pg.power for pg in mean_pgs], layout=(4,4),
    yticks=false, legend=false, size=(1600, 1200))
savefig(p, "graphics\\mean_periodograms.png")

p = plot(mean_pgs[1].freq/1e3, [pg.power for pg in mean_pgs], layout=(4,4),
    yticks=false, legend=false, xlims=(20, 40), xticks=20:40,
    gridalpha=0.3, size=(1600, 1200))
savefig(p, "graphics\\mean_periodograms_zoomed.png")

concatenated_clicks = [vec(waveforms[sub_i, :][labels .== i, :]') for i in 0:nlabels-1]
concatenated_sgs = [spectrogram(x, 512, fs=1/dt) for x in concatenated_clicks]
sg = concatenated_sgs[15]
heatmap(sg.time, sg.freq/1e3, 10log10.(sg.power), clim=(-100, -40),
    xlabel="Time (s)", ylabel="Freq (kHz)")
tiles4 = [heatmap(sg.time, sg.freq/1e3, 10log10.(sg.power),
    ylims=(0, 90),  clim=(-100, -40), xlabel="Time (s)", ylabel="Freq (kHz)", legend=false)
    for sg in concatenated_sgs]
p = plot(tiles4..., layout=(4, 4), size=(1600, 1200))
savefig(p, "graphics\\concatenated_spectrograms.png")

###############################################################################
# fitting decision tree
###############################################################################
using DecisionTree

p_train = 0.6
ii_train = StatsBase.sample(1:nsub, round(Int, p_train*nsub), replace=false)
ii_test = setdiff(1:nsub, ii_train)

model = build_forest( labels[ii_train], X[sub_i, :][ii_train, :])
nfoldCV_forest(labels[ii_train], X[sub_i, :][ii_train, :], 3)
save("output\\uninformed_classifier.jld2", Dict("model"=>model))

apply_forest(model, Vector{Float64}(features[1014, feature_names]))
plot(waveforms[1014, :])


###############################################################################
# choosing files for manual analysis
###############################################################################
using DataFramesMeta
using Dates

features[!, :label] = [apply_forest(model, Vector{Float64}(features[i, feature_names]))
    for i in 1:nrow(features)]

function clickfile2wavfile(filename)
    wavfile = split(basename(filename), "_clicks")[1]
    wavroot = "\\\\atlas.shore.mbari.org\\PAM_Archive_"
    y = wavfile[6:9]
    m = wavfile[10:11]
    return "$wavroot$y\\$m\\$wavfile"
end


label_files = @linq features |>
    where(:clicktime .< DateTime(2019, 1, 1, 0)) |>
    where(:label .!= 3) |> # exclude 38 kHz clicks
    by([:label, :filename], N = length(:clicktime)) |>
    transform(fileid = lpad.(1:length(:filename), 6, "0")) |>
    transform(filepath = clickfile2wavfile.(:filename)) |>
    transform(filename = basename.(:filepath))

diversity = @linq label_files |>
    by(:filename, shannon = -sum(:N ./ sum(:N) .* log.(:N ./ sum(:N))), totalN=sum(:N))
label_files = @linq label_files |>
    join(diversity, kind=:left, on=:filename) |>
    transform(P = :N ./ :totalN) |>
    where(:N .> 20, :P .> 0.2) # large sample of clicks, and large proportion of total

sort!(label_files, [:label, :N], rev=(false, true))

label_files = @linq vcat([g[2:2, :] for g in groupby(label_files, :label)]...) |>
    select(:label, :shannon, :N, :fileid, :filepath, :filename)

CSV.write("manual\\initial_files_for_tetyana.csv", label_files)

ftp_dir = "\\\\atlas.shore.mbari.org\\FTP\\pub\\urmy\\MARS_files"

for row in eachrow(label_files)
    outfile = "MARS_" * row.fileid * ".wav"
    cp(row.filepath, joinpath(ftp_dir, outfile))
    println(row.filepath, " >> ", joinpath(ftp_dir, outfile))
end
