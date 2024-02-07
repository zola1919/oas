import sys
import os


"""
A1-Part-2: Basic operations with audio

Write a function that reads an audio file and returns the minimum and the maximum values of the audio 
samples in that file. 

The input to the function is the wav file name (including the path) and the output should be two floating 
point values returned as a tuple.

If you run your code using oboe-A4.wav as the input, the function should return the following output:  
(-0.83486432, 0.56501967)"""
import numpy as np
from scipy.io import wavfile

def minMaxAudio(inputFile):

    fs, audio_data = wavfile.read(inputFile)
    audio_data=np.divide(audio_data,float(2**15-1),dtype=np.float32)
    min_val = np.min(audio_data)
    max_val = np.max(audio_data)

    return (min_val, max_val)


inputFile = 'oboe-A4.wav' 
min_max_values = minMaxAudio(inputFile)
print(min_max_values)