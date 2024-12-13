using DataFrames
using JLD2, FileIO
using CSV
using Dates

println("Listing raw .wav files...")
years = string.(2020:year(Dates.now()))
#months = lpad.(05:12, 2, "0")
month = "05"
#wav_root = "\\\\atlas.shore.mbari.org\\PAM_Archive"
click_root = "X:\\OdontoceteClicks\\clicks\\" .* years[1] .*"\\" .* month
wav_root = "\\\\atlas.shore.mbari.org\\PAM_Archive\\" .* years[1] .*"\\" .* month
#click_root = "C:\\clicks"


#years = string.(2016:2016)
#months = lpad.(10:10, 2, "0")
#day = "02"
#yeardirs = wav_root .* "\\" .* years .* "\\"
#raw_dirs = vcat([bd .* months .* "\\" for bd in yeardirs]...)
#raw_dirs = filter(isdir, raw_dirs)

wavfiles_vec = vcat([joinpath.(d, readdir(d)) for d in raw_dirs]...)
wavfiles_vec = readdir(wav_root)
wavfiles_vec = filter(f -> endswith(f, ".wav"), wavfiles_vec)
# Some duplicate wavfiles may exist in the wrong directories.  Eliminate them
# by making sure the date in the filename matches year/month directory
#wavfiles_vec = filter(f -> f[50:55] == f[37:40] * f[42:43], wavfiles_vec)
#wavfiles_vec = filter(f -> f[56:57] == day, wavfiles_vec)

# Next two lines calculate total raw data volume. Reading all the filesizes
# takes about 6 minutes, so they're commented out to avoid haning up the rest
# of the script.
# @time datasize = mapreduce(f -> filesize(f)/1e12, +, wavfiles_vec)
# println("Total raw data size: $(round(datasize)) TB")

wav_files = DataFrame(file=wavfiles_vec,
                       processed=false,
                       processed_time = DateTime(2020, 4, 1, 12, 00))
                       #processed_time = DateTime(2018, 5, 8, 12, 00))
allowmissing!(wav_files, :processed_time)
wav_files[:, :processed_time] .= missing

println("Listing detected click files...")
clickfiles_vec = []
for (root, dirs, files) in walkdir(click_root)
    ff = filter(f -> endswith(f, ".jld2"), files)
    append!(clickfiles_vec, joinpath.(click_root, root, ff))
end
processed_times = unix2datetime.(mtime.(clickfiles_vec))

println("Updating list of processed .wav files...")
fmt = "yyyymmdd_HHMMSS"
dts = [DateTime(basename(f)[6:20], fmt) for f in clickfiles_vec]
click_files = DataFrame(file = clickfiles_vec,
                        time = dts,
                        date = Date.(dts),
                        processed_time = processed_times)

# Find the .wav files which have corresponding detected-click files
a = [basename(f)[1:24] for f in clickfiles_vec]
b = [basename(f)[1:24] for f in wav_files.file]
ii = indexin(a, b)
wav_files[ii, :processed] .= true
wav_files[ii, :processed_time] .= processed_times

println("Saving file lists...")
CSV.write("output\\wav_files.csv", wav_files)
CSV.write("output\\click_files.csv", click_files)

println("Done!")
