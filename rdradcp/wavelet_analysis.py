import wavelets.kCwt


def plot_wavelet_spectrogram():

        kwavelet = wavelets.kCwt.kCwt(ppath, file, tunits)
        dj = 0.05  # Four sub-octaves per octaves
        s0 = -1  # 2 * dt                      # Starting scale, here 6 months
        J = -1  # 7 / dj                      # Seven powers of two with dj sub-octaves
        alpha = 0.5  # Lag-1 autocorrelation for white noise

        [wave, scales, freq, coi, fft, fftfreqs, iwave, power, fft_power, amplitude, phase] = \
            kwavelet.doSpectralAnalysis(title, "morlet", slevel, avg1, avg2, dj, s0, J, alpha)
        if debug:
            print("fftfreq=", fftfreqs)
            print("amplit=", amplitude)
            print("phase=", phase)


        ylabel_ts = "amplitude"
        yunits_ts = 'm'
        xlabel_sc = ""
        ylabel_sc = 'Period (%s)' % kwavelet.tunits
        # ylabel_sc = 'Freq (Hz)'
        sc_type = "period"
        # sc_type = "freq"
        # x_type = 'date'
        x_type = 'dayofyear'
        kwavelet.plotSpectrogram(ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2)
        ylabel_sc = 'Frequency ($s^{-1}$)'
        kwavelet.plotAmplitudeSpectrogram(ylabel_ts, yunits_ts, xlabel_sc, ylabel_sc, sc_type, x_type, val1, val2)

