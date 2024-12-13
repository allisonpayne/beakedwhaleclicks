
using Distributed

addprocs(2)
@everywhere using Pkg
@everywhere Pkg.activate(".")

@everywhere begin
    include("ClickDetection.jl")
    using FileIO
    using JLD2
    using DataFrames
    using CSV
    using Dates
    using DSP
    using Statistics
    using .ClickDetection
    using SharedArrays
end

@everywhere begin

    #click_root = "X:\\OdontoceteClicks\\clicks"
    click_root = "/Volumes/My Book/sandbox/"

    function peak_lead_ratio(click::ClickPointer, nsamples_lead)
        x = samples(click)
        t = times(click)
        x² = abs2.(x)
        return mean(x²[nsamples_lead:end]) / mean(x²[1:nsamples_lead])
    end

    function output_dir(wavfile)
        pathparts = split(wavfile, "\\")
        return joinpath(click_root, pathparts[end-2][end-3:end], pathparts[end-1])
    end

    function process_file(wavfile; noise_pct=0.6, noise_factor=200, span=2e-4, lead_snr=5)
        println("Reading $wavfile...")
        if filesize(wavfile) == 0
            println("Empty .wav file: $wavfile")
            return missing
        end

        audio = loadstreaming(wavfile)
        x = read(audio)
        close(audio)

        # preprocess raw audio before click detection
        x = convert.(Float32, x)
        fs = samplerate(x)
        butter = digitalfilter(Highpass(2e3, fs=fs), Butterworth(3))
        x = filtfilt(butter, x)

        # Notch filter to remove 50 kHz hum
        notch = iirnotch(50e3, 1e3, fs=fs)
        x = filt(notch, x)

        # Actual click detection
        x_teager = teager(x)
        noise_floor = quantile(vec(x_teager), noise_pct)
        thresh = noise_floor * noise_factor
        clicks = detect_clicks(x, fs, thresh, span)

        # eliminate "clicks" where high-teager detection was just a random spike
        nsamples_lead = floor(Int, span * fs)
        clicks = filter(c -> peak_lead_ratio(c, nsamples_lead) > lead_snr, clicks)

        # convert to DataFrame and save
        start_time = DateTime(basename(wavfile)[6:20], "yyyymmdd_HHMMSS")
        click_df = as_dataframe(clicks, start_time, wavfile)

        outdir = output_dir(wavfile)
        isdir(outdir) || mkpath(outdir)
        outfile = joinpath(outdir, basename(wavfile) * "_clicks.jld2")
        save(outfile, Dict("clicks" => click_df))

        return click_df
    end

end

function process_nfiles(wav_files, nfiles)
    rows_to_process = findall(.! wav_files.processed)[1:nfiles]
    pmap(process_file, wav_files.file[rows_to_process])
end

function process_hours(wav_files, hours)
    println("Processing for $(round(hours, digits=3)) hours")
    seconds = hours * 3600
    rows_to_process = findall(.! wav_files.processed)
    nfiles = length(rows_to_process)
    jobs = RemoteChannel(() -> Channel{Int}(nfiles))
    @async for i in rows_to_process
        put!(jobs, i)
    end
    completed = RemoteChannel(() -> Channel{Int}(nfiles))
    t0 = now()
    ncompleted = SharedArray{Int}(1)

    @everywhere function do_work(t0, seconds, all_wavfiles, jobs, completed, ncompleted)
        while (now() - t0).value / 1e3 < seconds
            i = take!(jobs)
            process_file(all_wavfiles[i])
            put!(completed, i)
            ncompleted[1] += 1
        end
    end
    @sync begin
        for p in workers()
            remote_do(do_work, p, t0, seconds, wav_files.file, jobs, completed, ncompleted)
        end
    end
    finalize(jobs)
    finalize(completed)
    sleep(seconds+1); println("Processed $(ncompleted[1]) files.")
    finalize(ncompleted)
end


function process_all(wav_files)
    nfiles = sum(.! wav_files.processed)
    process_nfiles(wav_files, nfiles)
end


# ii = findall(x -> x[50:57] ∈ ["20190201", "20190211"], wav_files[:file])
# wav_files[:processed] .= true
# wav_files[:processed][ii] .= false

include("list_files2.jl")
# wav_files = DataFrame(CSV.read("output/wav_files.csv"))
wav_files = DataFrame(CSV.read("wav_files.csv"))

# click_df = process_file(wav_files.file[end])
# load("F:\\clicks\\2020\\01\\MARS_20200101_063851.wav_clicks.jld2", "clicks")
#process_hours(wav_files, 24)
#process_nfiles(wav_files,4)
process_all(wav_files)

# n1 = sum(wav_files.processed)
# sum(.! wav_files.processed)
# @time process_hours(wav_files, 14)
# sum(wav_files.processed) - n1
