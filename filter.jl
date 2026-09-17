# Scratch comparison of filter padding strategies (odd-symmetric vs zero
# padding vs none). Not part of the pipeline; run it in a REPL and look at `fig`.
using CairoMakie   # in both shipped environments; swap for GLMakie if you want to zoom
using DSP

# Original signal: a simple oscillatory trend
t_min, t_max, dt = -40e3, -35e3, .4
t = t_min:dt:t_max  # Time vector

# Define periods in kyr
periods = [405, 100, 41, 21]

# Define amplitudes (arbitrary scaling)
amplitudes = [.3, 0.1, 0.5, 1]

# Generate signal as a sum of sinusoids
signal = sum(a * sin.(2π * t ./ p) for (a, p) in zip(amplitudes, periods))

# Add Gaussian noise
noise_level = 0.2
signal .+= noise_level * randn(length(t))

pad_len = length(signal)  # Length of extension

x = signal
# Sign-flipped mirror extension. NB: this is *not* what DSP.filtfilt does
# internally (that reflects about the end values, 2x[1] .- x[pad+1:-1:2]); it
# is just one more padding variant to compare against.
odd_sym_ext = vcat(-reverse(x[1:pad_len]), x, -reverse(x[end-pad_len+1:end]))

# Zero-padding extension
zero_pad_ext = vcat(zeros(pad_len), x, zeros(pad_len))

fs = 1/dt
nyq = fs/2
flow = (1/450) / nyq
fhigh = (1/400) / nyq
filter_order=4
# filter_order=2

bpass = digitalfilter(Bandpass(flow,fhigh),
                      Butterworth(filter_order))

simple_filter=filt(bpass,signal)
filtered_signal=filtfilt(bpass,signal)
filtered_signal_odd=filtfilt(bpass,odd_sym_ext)
filtered_signal_zero=filtfilt(bpass,zero_pad_ext)

bpass_fir = digitalfilter(Bandpass(flow, fhigh), FIRWindow(hamming(501)))
filtered_signal_hw = filtfilt(bpass_fir, signal)

bpass_fir = digitalfilter(Bandpass(flow, fhigh), FIRWindow(blackman(501)))
filtered_signal_bm = filtfilt(bpass_fir, signal)

bpass_fir = digitalfilter(Bandpass(flow, fhigh), FIRWindow(kaiser(501, 1.0)))
filtered_signal_kai = filtfilt(bpass_fir, signal)

# Create plots
fig = Figure()
ax = Axis(fig[1, 1])
# ax2 = Axis(fig[1, 2], title="Zero-Padding")

# Plot original signal and extensions
# lines!(ax, 1:length(odd_sym_ext), odd_sym_ext, label="Odd-Symmetric Extended")
lines!(ax, t, signal, label="Original Signal", alpha = 0.1)
lines!(ax, t, filtered_signal, label="Butterworth")
lines!(ax, t, simple_filter, label="filt Butter")
# lines!(ax, t, filtered_signal_hw, label="Hemmingworth")
# lines!(ax, t, filtered_signal_bm, label="Blackman")
# lines!(ax, t, filtered_signal_kai, label="Kaiser")
nb = length(signal)
# the padded series are [pad; signal; pad], so the original samples are nb+1:2nb
lines!(ax, t, filtered_signal_odd[nb+1:2nb], label="Odd-Symmetric Extended")
lines!(ax, t, filtered_signal_zero[nb+1:2nb], label="Zero-Padded Extended")
axislegend(ax)

fig
