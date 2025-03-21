#!/usr/bin/env python3
# Copyright (c) Max Planck Institute of Animal Behavior
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
"""
This script applies a quick & dirty event detection to audio files.
Implementation is loosely based on the one used here:
Denton, T., Wisdom, S., & Hershey, J. R. (2022).
Improving Bird Classification with Unsupervised Sound Separation.
2022 IEEE ICASSP, 636–640.
"""

import os
import sys
import torch
import shutil
import argparse
import torchaudio
import librosa
import numpy as np
import datetime as dt
from pathlib import Path
from scipy import signal
import soundfile as sf
from tqdm import tqdm

import matplotlib

import multiprocessing as mp

matplotlib.use("agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

sns.set_style("whitegrid")
butterworth_q = str(1 / np.sqrt(2))
custom_eq = [
    ["equalizer", "20", butterworth_q, "-20"],
    ["equalizer", "25", butterworth_q, "-20"],
    ["equalizer", "31", butterworth_q, "-20"],
    ["equalizer", "40", butterworth_q, "-20"],
    ["equalizer", "50", butterworth_q, "-20"],
    ["equalizer", "63", butterworth_q, "-20"],
    ["equalizer", "80", butterworth_q, "-9"],
    ["equalizer", "100", butterworth_q, "-8"],
    ["equalizer", "125", butterworth_q, "-5"],
    ["equalizer", "160", butterworth_q, "-3"],
    ["equalizer", "200", butterworth_q, "-2"]
]

# idea: use freqs above 1.5kHz - 2kHz

def hz_to_mel(frequencies, *, htk=False):
    """Convert Hz to Mels

    frequencies : number or np.ndarray [shape=(n,)] , float
        scalar or array of frequencies
    htk : bool
        use HTK formula instead of Slaney
    Returns
    -------
    mels : number or np.ndarray [shape=(n,)]
        input frequencies in Mels
    See Also
    --------
    mel_to_hz
    """

    frequencies = np.asanyarray(frequencies)

    if htk:
        return 2595.0 * np.log10(1.0 + frequencies / 700.0)

    # Fill in the linear part
    f_min = 0.0
    f_sp = 200.0 / 3

    mels = (frequencies - f_min) / f_sp

    # Fill in the log-scale part

    min_log_hz = 1000.0  # beginning of log region (Hz)
    min_log_mel = (min_log_hz - f_min) / f_sp  # same (Mels)
    logstep = np.log(6.4) / 27.0  # step size for log region

    if frequencies.ndim:
        # If we have array data, vectorize
        log_t = frequencies >= min_log_hz
        mels[log_t] = min_log_mel + np.log(frequencies[log_t] / min_log_hz) / logstep
    elif frequencies >= min_log_hz:
        # If we have scalar data, heck directly
        mels = min_log_mel + np.log(frequencies / min_log_hz) / logstep

    return mels


def mel_to_hz(mels, *, htk=False):
    """Convert mel bin numbers to frequencies

    mels : np.ndarray [shape=(n,)], float
        mel bins to convert
    htk : bool
        use HTK formula instead of Slaney
    Returns
    -------
    frequencies : np.ndarray [shape=(n,)]
        input mels in Hz
    See Also
    --------
    hz_to_mel
    """

    mels = np.asanyarray(mels)

    if htk:
        return 700.0 * (10.0 ** (mels / 2595.0) - 1.0)

    # Fill in the linear scale
    f_min = 0.0
    f_sp = 200.0 / 3
    freqs = f_min + f_sp * mels

    # And now the nonlinear scale
    min_log_hz = 1000.0  # beginning of log region (Hz)
    min_log_mel = (min_log_hz - f_min) / f_sp  # same (Mels)
    logstep = np.log(6.4) / 27.0  # step size for log region

    if mels.ndim:
        # If we have vector data, vectorize
        log_t = mels >= min_log_mel
        freqs[log_t] = min_log_hz * np.exp(logstep * (mels[log_t] - min_log_mel))
    elif mels >= min_log_mel:
        # If we have scalar data, check directly
        freqs = min_log_hz * np.exp(logstep * (mels - min_log_mel))

    return freqs


def mel_frequencies(n_mels=128, *, fmin=0.0, fmax=11025.0, htk=False):
    """Compute an array of acoustic frequencies tuned to the mel scale.
    The mel scale is a quasi-logarithmic function of acoustic frequency
    designed such that perceptually similar pitch intervals (e.g. octaves)
    appear equal in width over the full hearing range.
    Because the definition of the mel scale is conditioned by a finite number
    of subjective psychoaoustical experiments, several implementations coexist
    in the audio signal processing literature [#]_. By default, librosa replicates
    the behavior of the well-established MATLAB Auditory Toolbox of Slaney [#]_.
    According to this default implementation,  the conversion from Hertz to mel is
    linear below 1 kHz and logarithmic above 1 kHz. Another available implementation
    replicates the Hidden Markov Toolkit [#]_ (HTK) according to the following formula::
        mel = 2595.0 * np.log10(1.0 + f / 700.0).
    The choice of implementation is determined by the ``htk`` keyword argument: setting
    ``htk=False`` leads to the Auditory toolbox implementation, whereas setting it ``htk=True``
    leads to the HTK implementation.
    .. [#] Umesh, S., Cohen, L., & Nelson, D. Fitting the mel scale.
        In Proc. International Conference on Acoustics, Speech, and Signal Processing
        (ICASSP), vol. 1, pp. 217-220, 1998.
    .. [#] Slaney, M. Auditory Toolbox: A MATLAB Toolbox for Auditory
        Modeling Work. Technical Report, version 2, Interval Research Corporation, 1998.
    .. [#] Young, S., Evermann, G., Gales, M., Hain, T., Kershaw, D., Liu, X.,
        Moore, G., Odell, J., Ollason, D., Povey, D., Valtchev, V., & Woodland, P.
        The HTK book, version 3.4. Cambridge University, March 2009.
    n_mels : int > 0 [scalar]
        Number of mel bins.
    fmin : float >= 0 [scalar]
        Minimum frequency (Hz).
    fmax : float >= 0 [scalar]
        Maximum frequency (Hz).
    htk : bool
        If True, use HTK formula to convert Hz to mel.
        Otherwise (False), use Slaney's Auditory Toolbox.
    Returns
    -------
    bin_frequencies : ndarray [shape=(n_mels,)]
        Vector of ``n_mels`` frequencies in Hz which are uniformly spaced on the Mel
        axis.
    """

    # 'Center freqs' of mel bands - uniformly spaced between limits
    min_mel = hz_to_mel(fmin, htk=htk)
    max_mel = hz_to_mel(fmax, htk=htk)

    mels = np.linspace(min_mel, max_mel, n_mels)

    return mel_to_hz(mels, htk=htk)


def do_plot(spectrogram, idx1, idx2, unsmooted_data,
            smoothed_data, pi, sample_rate, hop_length):
    plt.close("all")
    f, ax = plt.subplots(nrows=2, ncols=1,
                         figsize=(10, 5))

    # The spectrogram
    ax[0].imshow(spectrogram[:, idx1:idx2],
                 origin="lower")
    spectrogram_shape = spectrogram[:, idx1:idx2].shape
    ax[0].set_aspect(spectrogram_shape[1] / 5 / spectrogram_shape[0])
    ax[0].axis(False)

    # The unsmoothed signal
    ax[1].plot(unsmooted_data,
               color="tab:blue", alpha=.125)

    # The smoothed signal
    ax[1].plot(smoothed_data,
               lw=3, c="tab:blue")

    # The found peaks
    ax[1].scatter(pi,
                  smoothed_data[pi],
                  color="tab:orange", s=100, zorder=10)

    ax[1].set_xlim(left=idx1, right=idx2)

    ticks_loc = ax[1].get_xticks().tolist()
    ax[1].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    dates = [x * hop_length / sample_rate for x in ticks_loc]
    x_datetime = [dt.datetime.fromtimestamp(d - 3600).strftime("%H:%M:%S") for d in dates]

    ax[1].set_xticklabels(x_datetime, rotation=45)
    ax[1].set_yticklabels([])
    ax[1].set_xlabel("Time in '%H:%M:%S'")
    sns.despine(top=True, right=True, left=True, bottom=True,
                offset=25, trim=True, ax=ax[1])
    f.tight_layout()
    return f


def get_log_mel_energies(wav, sample_rate, n_mels=96, fft_window=3.1e-2,
                         fft_hop=3.875e-3, scale="slaney", mel_power=1,
                         apply_bandpass=True, apply_background_correction=True,
                         apply_clamping=True, return_spectrogram=True,
                         mean_window_size=None, low_cutoff=0, high_cutoff=4000):
    """
    :param wav: The numpy / pytorch array that holds the audio data
    :param sample_rate: Sample rate of the audio data
    :param n_mels: Number of mel filter-banks. Default is 96
    :param fft_window: The length of an audio chunk when applying stft in seconds. Default is 31ms
    :param fft_hop: The hop length in seconds. Default is 3.875 ms
    :param scale: The Mel scale to use: Either "slaney" or "htk"
    :param mel_power: Exponent for the magnitude spectrogram, (must be > 0) e.g., 1 for energy, 2 for power, etc.
    :param apply_bandpass: (Optional) Apply a bandpass before converting to Mel
    :param apply_background_correction: (Optional) Subtract the median of every Mel band from the Mel spectrum
    :param apply_clamping: (Optional) Cut off values smaller than zero
    :param return_spectrogram: (Optional) Should we return the Mel Spectrogram
    :param mean_window_size: (Optional) The window size for the mean in the background correction,
                             if None; Use global mean and std
    :param low_cutoff: (Optional) The lower boundary for energy integration
    :param high_cutoff: (Optional) The upper boundary for energy integration
    """

    mel_energy = librosa.feature.melspectrogram(y=wav,
                                                sr=sample_rate,
                                                n_fft=int(fft_window * sample_rate),
                                                hop_length=int(fft_hop * sample_rate),
                                                win_length=int(fft_window * sample_rate),
                                                n_mels=n_mels,
                                                fmax=int(sample_rate / 2),
                                                norm="slaney",
                                                htk=False
                                                )
    log_mel_energy = librosa.power_to_db(
        mel_energy,
        ref=np.max
    )

    if apply_bandpass:
        if isinstance(wav, np.ndarray):
            wav = torch.tensor(wav, dtype=torch.float32)
        if wav.dim() == 1:
            wav = wav.unsqueeze(0)
        wav, _ = torchaudio.sox_effects.apply_effects_tensor(wav,
                                                              sample_rate,
                                                              custom_eq)

    mels = mel_energy
    mel_energies = log_mel_energy
    mel_energies -= mel_energies.min()
    mel_energies /= (mel_energies.max() + 1e-8)  # Scale between 0 and 1
    if isinstance(mel_energies, np.ndarray):
        mel_energies = torch.tensor(mel_energies, dtype=torch.float32)
    if apply_background_correction:
        if mean_window_size is None:
            threshold_per_mel_bin = mel_energies.mean(-1)[0] + mel_energies.std(-1)[0]
            if isinstance(threshold_per_mel_bin, np.ndarray):
                threshold_per_mel_bin = torch.tensor(threshold_per_mel_bin, dtype=torch.float32)

            mel_energies = (mel_energies.squeeze().t() - threshold_per_mel_bin.squeeze()).t()
        else:
            mel_energies_shape = mel_energies.shape[-1]
            mean_per_mel_bin = mel_energies.mean(-1)[0]
            std_per_mel_bin = mel_energies.std(-1)[0].squeeze()
            mean_filter = torch.full((1, 1, mean_window_size), 1 / mean_window_size)
            # Pad left and right with global mean
            padded_size = (n_mels, 1, mean_window_size + mel_energies_shape)
            padded_input = mean_per_mel_bin.view(-1, 1, 1) * torch.full(padded_size, 1.)
            half_ws = int(mean_window_size / 2)
            padded_input[:, 0, half_ws:-half_ws] = mel_energies
            # Calculate windowed average
            threshold_per_mel_bin = torch.nn.functional.conv1d(padded_input.cuda(),
                                                               mean_filter.cuda(),
                                                               padding="valid",
                                                               stride=1)[:, :, :-1].cpu().squeeze()
            mel_energies = ((mel_energies.squeeze() - threshold_per_mel_bin).t() - std_per_mel_bin).t()

    if apply_clamping:
        mel_energies = torch.clamp(mel_energies, min=0)

    if return_spectrogram:
        mel_spectrogram = mel_energies.clone().squeeze().cpu().detach().numpy()
    else:
        mel_spectrogram = None

    mel_bands = mel_frequencies(n_mels=n_mels, fmax=sample_rate)
    low_cutoff = np.argmin(np.abs(mel_bands - low_cutoff))
    high_cutoff = np.argmin(np.abs(mel_bands - high_cutoff))

    mel_energies = mel_energies.squeeze()[low_cutoff:high_cutoff].sum(0)  # use only between low and high cutoff
    mel_energies = mel_energies.squeeze().cpu().detach().numpy()

    return mel_energies, mel_spectrogram


def peak_detection(data, peak_kwargs,
                   indices=None,
                   apply_clamping=True):
    """
    Return the smoothed signal and the indices of the found peaks

    :param data: The numpy / pytorch array that holds the log mel, or whatever, energies
    :param peak_kwargs: A dictionary that holds the params for scipy's find_peaks routine
    :param indices: (Optional) Only evaluate 'wav' at these indices. A tuple (from, to) is expected here
    :param apply_clamping: (Optional) If True, then cut off negative values
    """
    if indices is None:
        idx1, idx2 = 0, len(np.squeeze(data))
    else:
        idx1, idx2 = indices[0], indices[1]

    # Smooth the signal
    print("Smoothing the data", flush=True)
    smoothed_data = signal.savgol_filter(data[idx1:idx2], polyorder=2,
                                         window_length=peak_kwargs["width"] * 4)

    if apply_clamping:
        smoothed_data = torch.clamp(torch.tensor(smoothed_data), min=0).squeeze().numpy()

    # Find the peaks
    print("Finding peaks", flush=True)
    peak_kwargs.pop("width")
    peak_indices_sign, peak_properties = signal.find_peaks(
        smoothed_data, **peak_kwargs
    )

    return smoothed_data, peak_indices_sign, peak_properties


def check_dir(p, remove=False):
    if os.path.isdir(p) and remove:
        shutil.rmtree(p)
    if not os.path.isdir(p):
        os.makedirs(p)


def main():
    parser = argparse.ArgumentParser("Simple energy based peak-event detection")
    parser.add_argument("tracks", nargs='+', type=Path, default=[], help='Path to tracks')
    parser.add_argument("-out", "--out_dir",
                        default=os.path.curdir,
                        type=str,
                        help="The output directory. (Default is the current directory)")
    parser.add_argument("-sr", "--sample_rate",
                        default=44100,
                        type=int,
                        help="The output sample rate. (Default: 44100 Hz)")
    parser.add_argument("-low", "--low_cutoff",
                        default=200,
                        type=int,
                        help="Lower boundary for energy integration. (Default is 2000 Hz)")
    parser.add_argument("-high", "--high_cutoff",
                        default=4000,
                        type=int,
                        help="Upper boundary for energy integration. (Default is 4000 Hz)")
    parser.add_argument("-ws", "--window_size",
                        default=5,
                        type=int,
                        help="We use this window size in seconds before and after some "
                             "found event for plotting. (Default: 5 s)")
    parser.add_argument("-ew", "--expected_width",
                        default=0.05,
                        type=float,
                        help="The expected width for the signal smoother in seconds. (Default: 0.05 s)")
    parser.add_argument("-pr", "--prominence",
                        default=0.4,
                        type=float,
                        help="The prominence for the peak finder. "
                             "See https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html "
                             "and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.peak_prominences.html"
                             "(Default: 0.4)")
    parser.add_argument("-p", "--plot",
                        default=False,
                        type=bool,
                        help="Should we do a plot of the loudest region? (Default: False)")
    parser.add_argument("-bp", "--band_pass",
                        default=True,
                        type=bool,
                        help="Should we apply our bandpass? (Default: True)")
    parser.add_argument("-pa", "--plot_all",
                        default=False,
                        type=bool,
                        help="Should we do a plot of all segments where something was detected (Default: False)")
    args = parser.parse_args()



    for track in args.tracks:
        process(track, args)

def process(track, args):
    if not track.exists():
        raise NameError(f"File {track} does not exist. If the path contains spaces, "
            "please try again after surrounding the entire path with quotes \"\".")
    fnamesub = os.path.basename(track)
    print(f"Working on: {track}")
    n_mels = 98  # --> number of mel bins
    fft_win = 0.032  # --> length of audio chunk when applying stft in seconds
    fft_hop = fft_win / 8  # --> hop_length in seconds
    hop_length = int(fft_hop * args.sample_rate)
    waveform, sample_rate = sf.read(track)
    out_base_dir = os.path.join(args.out_dir, os.path.basename(track)[:-4])
    check_dir(out_base_dir, remove=True)

    if sample_rate != args.sample_rate:
        # resample_args = (sample_rate, args.sample_rate, "kaiser_window",
        #                  64, 0.9475937167399596, 14.769656459379492)
        # resampler = torchaudio.transforms.Resample(  # Same as librosa "Kaiser Best", but faster
        #     *resample_args,
        #     dtype=waveform.dtype,
        # )
        # waveform = resampler(waveform)
        waveform = librosa.resample(waveform.squeeze(), orig_sr=sample_rate,
                                    target_sr=args.sample_rate, res_type="kaiser_best")

    print(f"generating log mel energies for {fnamesub}")
    log_mel_energies, mel_spectrogram = get_log_mel_energies(
        wav=waveform,
        sample_rate=args.sample_rate,
        n_mels=n_mels,
        fft_window=fft_win,
        fft_hop=fft_hop,
        scale="slaney",
        mel_power=1,
        apply_bandpass=args.band_pass,
        apply_background_correction=True,
        apply_clamping=True,
        return_spectrogram=True,
        mean_window_size=None,
        low_cutoff=args.low_cutoff,
        high_cutoff=args.high_cutoff
    )

    # Expected width
    expected_width = np.round(args.expected_width * args.sample_rate / hop_length).astype(int)
    peak_args = {
        "width": expected_width,
        "prominence": args.prominence
    }
    print(f"peak detection for {fnamesub}")
    smoothed_log_mel_energies, peak_indices, peak_properties = peak_detection(
        data=log_mel_energies,
        peak_kwargs=peak_args,
        indices=None,
        apply_clamping=True)

    out_dir = args.out_dir
    check_dir(out_dir)
    out_filename = os.path.join(out_dir, "{}.csv".format(
        os.path.basename(track)[:-len(".wav")]))
    print("Found {} events. Writing to: {}".format(len(peak_indices), out_filename))
    with open(out_filename, "w") as csv_file:
        csv_file.write("TimeInSeconds,Amplitude,WidthInSeconds\n")
        for pi, pp, prb, plb in tqdm(zip(peak_indices,
                                         peak_properties["prominences"],
                                         peak_properties["right_bases"],
                                         peak_properties["left_bases"]
                                         )):
            csv_file.write("{},{},{}\n".format(
                pi * hop_length / args.sample_rate,
                pp,
                (prb - plb) * hop_length / args.sample_rate
            )
            )

    if args.plot:
        print("Building a plot for the loudest found event +- {:02.02f}s.".format(args.window_size), flush=True)
        # We plot the region (+- n s) with the highest log mel energy (the loudest signal)
        loudest_idx = np.argmax(smoothed_log_mel_energies)
        n_seconds_in_log_mel_frames = int(args.window_size * args.sample_rate / hop_length)
        n_sec_before = max(  # Should not be negative
            0,
            loudest_idx - n_seconds_in_log_mel_frames
        )
        n_sec_after = min(  # Should not be larger than the file
            len(smoothed_log_mel_energies.squeeze()),
            loudest_idx + n_seconds_in_log_mel_frames
        )
        fig = do_plot(mel_spectrogram, n_sec_before, n_sec_after, log_mel_energies,
                      smoothed_log_mel_energies, peak_indices, args.sample_rate, hop_length)
        out_dir = os.path.join(out_base_dir, "images", "highest_energy")
        check_dir(out_dir)
        fig.savefig(os.path.join(out_dir, "{}_highest_energy_at_{:02.02f}s.png".format(
            os.path.basename(track),
            loudest_idx / args.sample_rate * hop_length)),
                    dpi=300)
        plt.close(fig)

    if args.plot_all:
        print("Building a plot for all events +- {:02.02f}s.".format(args.window_size), flush=True)

        all_windows = []
        n_seconds_in_log_mel_frames = int(args.window_size * args.sample_rate / hop_length)
        for pi in peak_indices:  # iterating over all found events
            n_sec_before = max(  # Should not be negative
                0,
                pi - n_seconds_in_log_mel_frames
            )
            n_sec_after = min(  # Should not be larger than the file
                len(smoothed_log_mel_energies.squeeze()),
                pi + n_seconds_in_log_mel_frames
            )
            all_windows.extend(np.arange(n_sec_before, n_sec_after))
        all_windows = list(set(all_windows))
        boundaries_for_windows = [0] + np.argwhere(
            np.diff(all_windows) > 1).squeeze().tolist() + [len(all_windows) - 1]

        for ba, bb in zip(boundaries_for_windows[:-1], boundaries_for_windows[1:]):
            # sanity check that all segments have no abrupt jumps
            sanity_one = np.diff(all_windows[ba + 1:bb + 1]).sum()
            sanity_two = np.ones(bb - ba).sum() - 1
            if sanity_one != sanity_two:
                print("Sanity check is off by {} for idx {} to {}".format(
                    np.abs(sanity_two - sanity_one), all_windows[ba], all_windows[bb]))
            segments_temp = np.arange(all_windows[ba], all_windows[bb])
            for i in range(np.ceil(len(segments_temp) / args.window_size).astype(int)):
                idx1 = int(i * 2 * n_seconds_in_log_mel_frames)
                idx2 = int(2 * (i + 1) * n_seconds_in_log_mel_frames)
                if idx1 < len(segments_temp):
                    idx2 = idx2 if idx2 < len(segments_temp) else len(segments_temp) - 1
                    # Check if some peak is within the current window
                    if any([segments_temp[idx1] < xx < segments_temp[idx2] for xx in peak_indices]):
                        fig = do_plot(mel_spectrogram, segments_temp[idx1], segments_temp[idx2], log_mel_energies,
                                      smoothed_log_mel_energies, peak_indices, args.sample_rate, hop_length)
                        out_dir = os.path.join(out_base_dir, "images", "all_events")
                        check_dir(out_dir)
                        fig.savefig(os.path.join(out_dir,
                                                 "{}_found_events_from_{:05.0f}s_to_{:05.0f}s.png".format(
                                                     os.path.basename(track),
                                                     np.round(segments_temp[idx1] / args.sample_rate * hop_length),
                                                     np.round(segments_temp[idx2] / args.sample_rate * hop_length))
                                                 ),
                                    dpi=300)
                        plt.close(fig)


if __name__ == "__main__":
    main()
