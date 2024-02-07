import os
import sys
import numpy as np
from scipy.signal import get_window
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import stft
import utilFunctions as UF
from scipy.io import wavfile
eps = np.finfo(float).eps
from scipy.io import wavfile
"""
A4-Part-3: Computing band-wise energy envelopes of a signal

Write a function that computes band-wise energy envelopes of a given audio signal by using the STFT.
Consider two frequency bands for this question, low and high. The low frequency band is the set of 
all the frequencies between 0 and 3000 Hz and the high frequency band is the set of all the 
frequencies between 3000 and 10000 Hz (excluding the boundary frequencies in both the cases). 
At a given frame, the value of the energy envelope of a band can be computed as the sum of squared 
values of all the frequency coefficients in that band. Compute the energy envelopes in decibels. 

Refer to "A4-STFT.pdf" document for further details on computing bandwise energy.

The input arguments to the function are the wav file name including the path (inputFile), window 
type (window), window length (M), FFT size (N) and hop size (H). The function should return a numpy 
array with two columns, where the first column is the energy envelope of the low frequency band and 
the second column is that of the high frequency band.

Use stft.stftAnal() to obtain the STFT magnitude spectrum for all the audio frames. Then compute two 
energy values for each frequency band specified. While calculating frequency bins for each frequency 
band, consider only the bins that are within the specified frequency range. For example, for the low 
frequency band consider only the bins with frequency > 0 Hz and < 3000 Hz (you can use np.where() to 
find those bin indexes). This way we also remove the DC offset in the signal in energy envelope 
computation. The frequency corresponding to the bin index k can be computed as k*fs/N, where fs is 
the sampling rate of the signal.

To get a better understanding of the energy envelope and its characteristics you can plot the envelopes 
together with the spectrogram of the signal. You can use matplotlib plotting library for this purpose. 
To visualize the spectrogram of a signal, a good option is to use colormesh. You can reuse the code in
sms-tools/lectures/4-STFT/plots-code/spectrogram.py. Either overlay the envelopes on the spectrogram 
or plot them in a different subplot. Make sure you use the same range of the x-axis for both the 
spectrogram and the energy envelopes.

NOTE: Running these test cases might take a few seconds depending on your hardware.

Test case 1: Use piano.wav file with window = 'blackman', M = 513, N = 1024 and H = 128 as input. 
The bin indexes of the low frequency band span from 1 to 69 (69 samples) and of the high frequency 
band span from 70 to 232 (163 samples). To numerically compare your output, use loadTestCases.py
script to obtain the expected output.

Test case 2: Use piano.wav file with window = 'blackman', M = 2047, N = 4096 and H = 128 as input. 
The bin indexes of the low frequency band span from 1 to 278 (278 samples) and of the high frequency 
band span from 279 to 928 (650 samples). To numerically compare your output, use loadTestCases.py
script to obtain the expected output.

Test case 3: Use sax-phrase-short.wav file with window = 'hamming', M = 513, N = 2048 and H = 256 as 
input. The bin indexes of the low frequency band span from 1 to 139 (139 samples) and of the high 
frequency band span from 140 to 464 (325 samples). To numerically compare your output, use 
loadTestCases.py script to obtain the expected output.

In addition to comparing results with the expected output, you can also plot your output for these 
test cases.You can clearly notice the sharp attacks and decay of the piano notes for test case 1 
(See figure in the accompanying pdf). You can compare this with the output from test case 2 that 
uses a larger window. You can infer the influence of window size on sharpness of the note attacks 
and discuss it on the forums.
"""
def computeEngEnv(inputFile, window, M, N, H):
    """
    Inputs:
            inputFile (string): input sound file (monophonic with sampling rate of 44100)
            window (string): analysis window type (choice of rectangular, triangular, hanning, 
                hamming, blackman, blackmanharris)
            M (integer): analysis window size (odd positive integer)
            N (integer): FFT size (power of 2, such that N > M)
            H (integer): hop size for the stft computation
    Output:
            The function should return a numpy array engEnv with shape Kx2, K = Number of frames
            containing energy envelop of the signal in decibles (dB) scale
            engEnv[:,0]: Energy envelope in band 0 < f < 3000 Hz (in dB)
            engEnv[:,1]: Energy envelope in band 3000 < f < 10000 Hz (in dB)
    """
    
    ### your code here
        # read audio file
    fs, x = wavfile.read(inputFile)
       
       # get respective window
    w = get_window(window, M)
       
       # call the sftf with parameters
    xmX, xpX = stft.stftAnal(x, w, N, H)
       
    cutoff_f = 3000
    ceil_f = 10000
       
       # get the bin indexes of the 3000 and 10000 frequencies
    cutoff_bin = int(np.ceil(cutoff_f * N / fs))
    ceil_bin =  int(np.ceil(ceil_f * N / fs))
       
       # initialize the low and high envelopes
    low_env = np.zeros(len(xmX))
    high_env = np.zeros(len(xmX))
       
    for i in np.arange(len(xmX)):
           # for each magnitude spectrum (one per Hop) we obtain the
           # respective bins for each band frequency
           low_band = xmX[i][1:cutoff_bin]
           high_band = xmX[i][cutoff_bin:ceil_bin]
           
           # the bands are in dB, so we use the inverse of the dB
           # transformation and over the linear value we compute
           # the energy of the band
           energy_low_computed = np.sum(abs(10**(low_band/20))**2)
           energy_high_computed = np.sum(abs(10**(high_band/20))**2)
           
           # transform it back to dB
           energy_low_db = 10*np.log10(energy_low_computed)
           energy_high_db = 10*np.log10(energy_high_computed)
           
           # save the energy value for that window in the respective
           # array to be returned
           low_env[i] = energy_low_db
           high_env[i] = energy_high_db
       
    return xmX, low_env, high_env
   



def plot_spectrogram(xmX, M, N, H, fs, low_env, high_env, window,title = "Spectrogram vs Energy Envelopes"):
    # part of the code is from sms-tools/lectures/04-STFT/plots-code/spectrogram.py:
    plt.figure(2, figsize=(15, 10))

    plt.subplot(211)
    numFrames = int(xmX[:,0].size)
    frmTime = H*np.arange(numFrames)/fs                            
    binFreq = np.arange(N/2+1)*float(fs)/N                         
    plt.pcolormesh(frmTime, binFreq, np.transpose(xmX))
    plt.title('xmX, M=' + str(M) + ', N=' + str(N) + ', H=' + str(H))
    plt.autoscale(tight=True)
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (s)')



    plt.subplot(212)
    
    # same X_axis as in the spectrogram
    numFrames = int(len(low_env))
    frmTime = H*np.arange(numFrames)/fs                            

    # plot the respective arrays
    plt.plot(frmTime, low_env)
    plt.plot(frmTime, high_env)

    # plot params
    plt.title(title)
    plt.ylabel('Amplitude (dB)')
    plt.xlabel('Time (s)')
    plt.legend(['Low Frequency','High Frequency'])
    plt.autoscale(tight=True)
    plt.tight_layout()
    plt.savefig('{}.png'.format(title+"_"+window))
    plt.show()
  
    

    
if __name__ == "__main__":   
    input_file = '../../../sounds/piano.wav'
    window = 'blackman'
    M = 513
    N = 1024
    H = 128
    
    fs = 44100
    
    xmX, low_env, high_env = computeEngEnv(input_file, window, M, N, H)
    
    plot_spectrogram(xmX, M, N, H, fs,low_env, high_env,window)
