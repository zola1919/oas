import os
import sys
import numpy as np
import math
from scipy.signal import get_window
import matplotlib.pyplot as plt
from scipy.io import wavfile
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../../software/models/'))

import stft
import utilFunctions as UF
eps = np.finfo(float).eps


"""
A4-Part-2: Measuring noise in the reconstructed signal using the STFT model 

Write a function that measures the amount of noise introduced during the analysis and synthesis of a 
signal using the STFT model. Use SNR (signal to noise ratio) in dB to quantify the amount of noise. 
Use the stft() function in stft.py to do an analysis followed by a synthesis of the input signal.

A brief description of the SNR computation can be found in the pdf document (A4-STFT.pdf, in Relevant 
Concepts section) in the assignment directory (A4). Use the time domain energy definition to compute
the SNR.

With the input signal and the obtained output, compute two different SNR values for the following cases:

1) SNR1: Over the entire length of the input and the output signals.
2) SNR2: For the segment of the signals left after discarding M samples from both the start and the 
end, where M is the analysis window length. Note that this computation is done after STFT analysis 
and synthesis.

The input arguments to the function are the wav file name including the path (inputFile), window 
type (window), window length (M), FFT size (N), and hop size (H). The function should return a python 
tuple of both the SNR values in decibels: (SNR1, SNR2). Both SNR1 and SNR2 are float values. 

Test case 1: If you run your code using piano.wav file with 'blackman' window, M = 513, N = 2048 and 
H = 128, the output SNR values should be around: (67.57748352378475, 304.68394693221649).

Test case 2: If you run your code using sax-phrase-short.wav file with 'hamming' window, M = 512, 
N = 1024 and H = 64, the output SNR values should be around: (89.510506656299285, 306.18696700251388).

Test case 3: If you run your code using rain.wav file with 'hann' window, M = 1024, N = 2048 and 
H = 128, the output SNR values should be around: (74.631476225366825, 304.26918192997738).

Due to precision differences on different machines/hardware, compared to the expected SNR values, your 
output values can differ by +/-10dB for SNR1 and +/-100dB for SNR2.
"""

def computeSNR(inputFile, window, M, N, H):
    # Za ulazni signal
    fs, x = wavfile.read(inputFile)

    y = stft.stft(x, w=get_window(window, M), N=N, H=H)


    # Compute the time-domain energy of the original signal
    original_energy = np.sum(x**2)

    # Compute the time-domain energy of the noise signal (original - reconstructed)
    noise_signal = x - y  # Provera da se duzine podudaraju
    noise_energy = np.sum(noise_signal**2)

    # Obrazac za racunanje SNR (Signal to noise ratio) u dB
    SNR1 = 10 * np.log10(original_energy / noise_energy)
    
    # Racunanje SNR2 za segment koji je ostao nakon odbacivanja M uzoraka i sa pocetka i sa kraja
    x_trimmed = x[M:-M]  #Odbacite tih M uzoraka
    y_trimmed = y[M:-M]
    original_energy_trimmed = np.sum(x_trimmed**2)
    noise_signal_trimmed = x_trimmed - y_trimmed
    noise_energy_trimmed = np.sum(noise_signal_trimmed**2)
    SNR2 = 10 * np.log10(original_energy_trimmed / noise_energy_trimmed)

    return (SNR1, SNR2)
if __name__ == "__main__": 
    inputFile = "../../../sounds/piano.wav"
    window = "blackman"
    M = 512
    N = 2048
    H = 128
    snr_values = computeSNR(inputFile, window, M, N, H)
    print("SNR1:", snr_values[0])
    print("SNR2:", snr_values[1])