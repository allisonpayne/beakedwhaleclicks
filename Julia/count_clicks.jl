using CSV
using JLD2
using UMAP
using ConcaveHull
using FileIO
using Arrow
using DSP
using Statistics
using Suppressor
using DataFrames, DataFramesMeta
using AxisArrays
using SampledSignals
using Dates
using Alert

logfile = joinpath(@__DIR__, "../data/log.txt")

include(joinpath(@__DIR__, "deimos_detector.jl"))
include(joinpath(@__DIR__, "distance_metrics.jl"))
include(joinpath(@__DIR__, "rissos_pwsd.jl"))

const fs = 256e3
n_neighbors = 6
min_dist = 0.001
umap_model = load(joinpath(@__DIR__, "../data/umap_model.jld2"), "umap_model")
cluster_hulls = load(joinpath(@__DIR__, "../data/cluster_hulls.jld2"), "cluster_hulls")

function pad(x, n=512)
    if length(x) > n
        return x[1:n]
    else
        return [x; zeros(n - length(x))]
    end
end

function read_waveforms(clickfile)
    clicks = nothing
    try
        clicks = @suppress load(clickfile, "clicks");
    catch
        msg = "!!!Failed to open $clickfile \n"
        open(logfile, "a") do f
            write(f, msg)
        end
        println(msg)
        return missing, missing
    end

    if :start_datetime in clicks.colindex.names
        i_time = clicks.colindex.lookup[:start_datetime]
    elseif :start in clicks.colindex.names
        i_time = clicks.colindex.lookup[:start]
    else # no data in file
        return missing, missing
    end
    
    i_wave = clicks.colindex.lookup[:waveform]
    # no clicks detected in this file
    if all(ismissing.(clicks.columns[i_wave]))
        return missing, missing
    end
    waveforms = Vector.(clicks.columns[i_wave])
    waveforms = [pad(w) for w in waveforms]
    waveforms = [(w .- mean(w)) ./ std(w) for w in waveforms]
    @assert all(var.(waveforms) .≈ 1)
    times = clicks.columns[i_time]
    return times, waveforms
end

function calculate_spectra(waveforms, n=256, fs=256e3)
    spectra = map(waveforms) do w
        x = [zeros(256); w; zeros(256)]
        welch_pgram(x, n, fs=fs)
    end
    return spectra
end

function embed(spectra)
    S = hcat([s.power for s in spectra]...)
    if length(spectra) <= n_neighbors
        k = div(n_neighbors, size(S, 2)) + 1
        S = repeat(S, outer=(1, k))
    end
    E = UMAP.transform(umap_model, S, metric=EarthMovers(), 
        n_neighbors=n_neighbors, min_dist=min_dist)
    E = E[:, 1:length(spectra)]
    return E
end

function classify(embedding, hulls)
    c = map(eachcol(embedding)) do em
        c = findfirst(h -> in_hull(em, h), hulls)
        isnothing(c) ? 0 : c
    end
    return c
end

function dolphin_labels(spectra, clusters)
    label = fill("NA", length(clusters))
    for i in eachindex(label)
        if clusters[i] == 1
            s = spectra[i]
            if rissos(s.freq, s.power)
                label[i] = "Risso's"
            elseif pwsd(s.freq, s.power)
                label[i] = "PWSD"
            else
                label[i] = "UID"
            end
        end
    end
    return label
end


function read_and_count(filename)
    println(filename)
    times, waveforms = read_waveforms(filename)
    if ismissing(waveforms)
        return nothing
    end
    spectra = calculate_spectra(waveforms)
    i_bio = findall(s -> !is_deimos(s.freq, s.power), spectra)
    if length(i_bio) == 0
        return nothing
    end
    spectra = spectra[i_bio]
    times = times[i_bio]
    embedding = embed(spectra)
    cluster = classify(embedding, cluster_hulls)
    dolphin = dolphin_labels(spectra, cluster)

    clickrate = DataFrame(
        clicktime = times,
        cluster = cluster,
        dolphin = dolphin
    )
    clickrate = @chain clickrate begin
        @transform(:datetime_round = round.(:clicktime, Minute(10)))
        @by([:datetime_round, :cluster, :dolphin], :nclicks = length(:clicktime))
        @orderby(:datetime_round, :cluster)
    end
    return clickrate
end


click_dir = "/media/sam/Sam01/clicks"
datadir = joinpath(@__DIR__, "..", "data")

# fn = "/media/sam/Sam01/clicks/2019/07/MARS_20190718_082900.wav_clicks.jld2"
# clicks = load(fn, "clicks");
# times, waveforms = read_waveforms(fn)
# ismissing(waveforms)
# spectra = calculate_spectra(waveforms)
# i_bio = findall(s -> !is_deimos(s.freq, s.power), spectra)
# spectra = spectra[i_bio]
# times = times[i_bio]
# embedding = embed(spectra)

# using StatsPlots
# p = scatter(embedding[1, :], embedding[2, :])
# for h in cluster_hulls
#     plot!(p, h)
# end
# p

# cluster = classify(embedding, cluster_hulls)
# clickrate = DataFrame(clicktime = times, cluster = cluster)
# clickrate = @chain clickrate begin
#     @transform(:datetime_round = round.(:clicktime, Minute(10)))
#     @by([:datetime_round, :cluster], :nclicks = length(:clicktime))
#     @orderby(:datetime_round, :cluster)
# end

# f = read_and_count(fn)

yearmonths = Date(2019, 01):Month(1):Date(2020, 01)

for date in yearmonths[11:end]
    yr = string(year(date))
    mnth = lpad(month(date), 2, "0")
    dir = joinpath(click_dir, yr, mnth)
    clickfiles = sort(readdir(dir, join=true))

    clickrate = map(clickfiles) do filename
        open(logfile, "a") do f
            write(f, filename * "\n")
        end
        read_and_count(filename)
    end
    
    clickrate = filter(!isnothing, clickrate)
    clickrate = reduce(vcat, clickrate)
    clickrate = @by(clickrate, [:datetime_round, :cluster, :dolphin], 
        :nclicks = sum(:nclicks))
    CSV.write(joinpath(datadir, "clickrate_umap_$(yr)$(mnth).csv"), clickrate)
end
alert("Click processing done!")
