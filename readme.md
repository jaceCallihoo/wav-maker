# Wave File Generator

## what it does
Generate audio data in the wave file format based on combined audio waves

## how it does it
Audio generation is handled through a streaming approach, where each buffer of samples is computed on demand and written to the output file or audio device before the next chunk is processed. This ensures low memory overhead and scalable audio production.

implemented multiple wave types
sin, square, sawtooth, triangle

imlemented linear transformation 
used runtime polymorphism to allow for both constant and linear transformations to aplitude and frequency

implemented linear transformatinos for frequency by calculating the average area under the curve for the given keyframes

accounted for memory endianness

## wave functions explained
each function accepts a time and frequency and returns the amlpitude of the wave at that time. 
The time is always relative to 0 as this prevents phase cancelation when adding waves together

### References
https://web.archive.org/web/20250610165025/http://soundfile.sapp.org/doc/WaveFormat/

