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


sound2 <- clip_by_sec(sound, fs, 0, 0.3)
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
