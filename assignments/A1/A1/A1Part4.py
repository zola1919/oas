import sys
import os

from A1Part3 import hopSamples
import numpy as np
from scipy.io import wavfile

"""
A1-Part-4: Downsampling audio: Changing the sampling rate

One of the required processes to represent a signal inside a computer is sampling. The sampling rate is the number of samples obtained in one second when sampling a continuous analog signal to a discrete digital signal. As mentioned earlier, most of the time we will be working with wav audio files that have a sampling rate of 44100 Hz, which is a typical value. For some applications, changing the sampling rate of an audio signal can be necessary. This optional part shows how to do this, from a higher sampling rate to a lower one.

Complete the function downsampleAudio(inputFile,M) in the file A1Part4.py so that given an audio file, it applies downsampling by a factor of M and create a wav audio file <input_name>_downsampled.wav at a lower sampling rate.

In Part1 you learned how to read a wav file and the function from Part3 can be used to perform the downsampling of a signal contained in an array. To create a wav audio file from an array, you can use the wavwrite function from the utilFunctions module. Be careful with the sampling rate parameter since it should be different from that of the original audio.

You can test your code using the file `vibraphone-C6.wav' and a downsampling factor of M=16. 
Listen to the `vibraphone-C6_downsampled.wav' sound. What happened to the signal?
How could we avoid damaging the signal when downsampling it?
You can find some related information in https://en.wikipedia.org/wiki/Decimation_%28signal_processing%29.
"""

def downsampleAudio(inputFile, M):
    """
    Inputs:
        inputFile: file name of the wav file (including path)
        	M: downsampling factor (positive integer)
    """
    fs, audio_data = wavfile.read(inputFile)
    print(audio_data.size/fs)
    fs = fs // M;
    y = hopSamples(audio_data, M)
    print(y.size/fs)
    name=inputFile.split("/")[-1][0:-4]+'_downsampled.wav'
    wavfile.write(name,fs,y)
    return(name)
    ## Your code here

inputFile = 'vibraphone-C6.wav'  # Add '.wav' to the filename
result = downsampleAudio(inputFile,5)
print(result)

