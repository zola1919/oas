import sys
import os

import numpy as np
from scipy.io import wavfile

def readAudio(inputFile):
    """
    Input:
    inputFile: the path to the wav file
    Output:
    The function should return a numpy array that
    contains 10 samples of the audio.
    """

    fs, audio_data = wavfile.read(inputFile)
    audio_data=np.divide(audio_data,float(2**15-1),dtype=np.float32)
    first_sample = 50000 
    last_sample = first_sample + 10
    ten_samples = audio_data[first_sample:last_sample]
    ten_samples=np.append(ten_samples,ten_samples.dtype)
    return list(ten_samples)

inputFile = 'piano.wav' 
result = readAudio(inputFile)
!play result
print(result)