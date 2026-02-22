# audio editor tool

- should create .wav files
- should edit .wav files? (maybe later)
- should be able to write to files
- how do I choose the endianness of a written byte



### References
https://web.archive.org/web/20250610165025/http://soundfile.sapp.org/doc/WaveFormat/


```
// ChunkID          constant
// ChunkSize        derived
// Format           constant
// Subchunk1ID      constant
// Subchunk1Size    constant maybe?
// AudioFormat      constant (1 for PCM)
// NumChannels      variable (1 = mono, 2 = stereo)
// SampleRate       variable (8000, 44100, etc)
// ByteRate         derived
// BlockAlign       derived
// BitsPerSample    variable (8 = 8 bits, 16 = 16 bits)
// ExtraParmasSize  other (if PCM, then doesn't exist)
// ExtraParams      other (space for extra parameters)
// Subchunk2ID      constant ("data")
// Subchunk2Size    derived (NumSamples * NumChannels * BitsPerSample/8)
// Data             variable
```


## Wave Math

Time (duration) x sample rate is used to determine the number of samples
For each sample, sample index x sample rate is used to determine time
Time x frequency used to determine amplitude



