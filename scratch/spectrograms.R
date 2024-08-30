library(tuneR)
library(tidyverse)
library(signal)
library(oce)

data <- readWave("data/MARS_20180613_142912.wav")
sound <- data@left
duration <- length(sound) / data@samp.rate
fs <- data@samp.rate

#plot waveform 

#how to decimate:
# sound2 <- sound[seq(1, length(sound), by = 1000)]
#fs2 <- fs / 1000

# clip a minute out
clip_by_sec <- function(x, f, start, end) {
  n_x <- as.numeric(length(x))
  f <- as.numeric(f)
  stopifnot(length(start) == 1,
            length(end) == 1,
            start >= 0,
            end >= 0,
            start <= length(x) * f,
            end <= length(x) * f)
  start_idx <- 1 + start * f
  end_idx <- 1 + end * f
  x[start_idx:end_idx]
}


sound2 <- clip_by_sec(sound, fs, 0, 300)
fs2 <- fs 

# Takes so looooong
#demean so the y axis is centered on zero
# sound2 <- sound2 - mean(sound2)
# plot(sound2, type = "l")

#create spectrogram
spec <- specgram(x = sound2, 
                 n = 1024, 
                 Fs = fs2, 
                 window = 256, 
                 overlap = 128)

#discard phase information
p <- abs(spec$S)
#normalize
p <- p / max(p)
#config time axis
t <- spec$t

#plot spec
# imagep(x = t, 
#        y = spec$f, 
#        z = t(p), 
#        col = oce.colorsViridis, 
#        ylab = "Frequency", 
#        xlab = "Time", 
#        drawPalette = TRUE, 
#        decimate = TRUE)

# zoom in on 20 - 120 kHz
flim <- c(20e3, 120e3)
plim <- range(p[between(spec$f, flim[1], flim[2]), ])
imagep(x = t, 
       y = spec$f, 
       z = t(p), 
       col = oce.colorsViridis, 
       ylab = "Frequency", 
       xlab = "Time", 
       #xlim = c(80, 90),
       ylim = flim,
       zlim = plim,
       drawPalette = TRUE, 
       decimate = FALSE)
clicks <- c(0.0341, 0.0776, 0.0901, 0.1446, 0.1889, 0.1942, 0.1987, 0.2514)
abline(v = clicks, col = "red", alpha = 0.2)


# Mean and center frequency of those clicks?
clicks_mat <- R.matlab::readMat("data/training/yearly/2018/clicks/MARS_20180613_142912.wav_clicks.mat")
# IS IT THE PEAK SUMMED ACROSS TIME??
peak_freq <- function(p, t, freq, t_lim, freq_lim) {
  in_t <- t_lim[1] <= t & t_lim[2] >= t
  in_freq <- freq_lim[1] <= freq & freq_lim[2] >= freq
  freq <- freq[in_freq]
  freq[which.max(rowSums(p[in_freq, in_t, drop = FALSE]))]
}
# IS ARITHMETIC MEAN EVEN RIGHT??
center_freq <- function(p, t, freq, t_lim, freq_lim) {
  in_t <- t_lim[1] <= t & t_lim[2] >= t
  in_freq <- freq_lim[1] <= freq & freq_lim[2] >= freq
  p <- p[in_freq, in_t, drop = FALSE]
  p_norm <- p / sum(p)
  freq <- freq[in_freq]
  freq_mat <- matrix(freq, 
                     nrow = length(freq), 
                     ncol = ncol(p_norm), 
                     byrow = TRUE)
  sum(p_norm * freq_mat)
}
p <- abs(spec$S)
p <- p / max(p)
clicks_df <- tibble(
  start = as.numeric(clicks_mat$start),
  stop = as.numeric(clicks_mat$stop),
  dur = stop - start
) %>% 
  dplyr::filter(stop < 300,
                dur > 0.001) %>% 
  mutate(
    peak_freq = map2_dbl(
      start, stop,
      \(t1, t2) peak_freq(p, spec$t, spec$f, 
                          t_lim = c(t1, t2),
                          freq_lim = c(5e3, 120e3))),
    center_freq =  map2_dbl(
      start, stop,
      \(t1, t2) center_freq(p, spec$t, spec$f, 
                            t_lim = c(t1, t2),
                            freq_lim = c(5e3, 120e3)))
  )

t <- spec$t
f <- spec$f
imagep(x = t, 
       y = f, 
       z = t(p), 
       col = oce.colorsViridis, 
       ylab = "Frequency", 
       xlab = "Time", 
       xlim = c(192.1, 192.3),
       ylim = flim,
       zlim = plim,
       drawPalette = TRUE, 
       decimate = FALSE)
# peak freq
points(x = mean(clicks_df$start[1], clicks_df$stop[1]),
       y = clicks_df[1, "peak_freq"],
       col = "red")
points(x = mean(clicks_df$start[1], clicks_df$stop[1]),
       y = clicks_df[1, "center_freq"],
       col = "orange")



# Signal frequency power --------------------------------------------------

# First: high pass filter above 10^3.5 Hz
nyquist <- fs / 2
W <- c(10^4, 10^5) / nyquist
bf <- butter(4, W, type = "pass")
sound_high <- filtfilt(bf, sound)

window_size_s <- 0.1
window_size_i <- window_size_s * fs
myclick_win <- as.integer(myclick_t * fs + c(-1, 1) * window_size_i / 2)
myclick_signal <- sound_high[myclick_win[1]:myclick_win[2]]
myclick_fft <- fft(myclick_signal)
n <- length(myclick_signal)
f <- (0:(n - 1)) * (fs / n)
power = abs(myclick_fft)^2 / n
in_range <- f >= 10^4 & f <= 10^5
plot(f[in_range], power[in_range], type = "l")

peak_freq <- f[which.max(power)]
center_freq <- weighted.mean(f, power)



