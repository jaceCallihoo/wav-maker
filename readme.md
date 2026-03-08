# Wave File Generator

A C++20 audio synthesizer that generates WAV files from programmatically defined waveforms. Audio is produced through a streaming pipeline — samples are computed on demand and written in fixed-size chunks, keeping memory usage constant regardless of output duration.

## Demos
```
Headphone Warning
```

### Constant Waves

#### Sin Wave 440hz
https://github.com/user-attachments/assets/9b040c27-f6b2-4277-b0f9-d5acd06efa09

#### Square Wave 440hz
https://github.com/user-attachments/assets/1170990d-4e6e-4799-9a3a-b60f9bcf8267

### Triangle Wave 440hz
https://github.com/user-attachments/assets/f66e8c35-5ffd-4910-ba98-0b7b540222d8

#### Sawtooth Wave 440hz
https://github.com/user-attachments/assets/6aa29b20-d7c6-474b-95e9-c7f5112e4c4d

### Linear Amplitude Transforms
#### Sin Wave Ease In 110hz
https://github.com/user-attachments/assets/c7d67b9a-8851-4240-9508-7195f8c98422

#### Triangle Wave Ease In Out 220hz
https://github.com/user-attachments/assets/187174e0-7f07-4485-8450-66090c6a0713

### Lineare Frequency Transforms
#### Square Wave 440hz to 880hz
https://github.com/user-attachments/assets/6c562edc-668f-443c-867c-8b151b072e42

### Combined
### Sin Wave x3, Ease In Out, 110hz to 220hz, 220hz to 440hz, 440hz to 880hz 
https://github.com/user-attachments/assets/0ea03731-be17-4ef1-9b68-67cbb083f5c2


## Features

- **Four waveform types** — sine, square, sawtooth, and triangle
- **Keyframe-based modulation** — smoothly vary amplitude and frequency over time using percent-based keyframes
- **Correct frequency integration** — pitch sweeps are calculated by integrating the area under the frequency curve (trapezoid rule), not just interpolating the value, producing accurate pitch output
- **Runtime-polymorphic transformations** — `Linear` and `Constant` transformation types share an abstract interface, allowing them to be composed and swapped freely
- **Sound mixing** — multiple sounds with independent delays and durations are summed and averaged to prevent clipping
- **WAV format output** — generates a spec-compliant 44-byte header followed by 16-bit little-endian PCM samples
- **Streaming architecture** — audio is generated and written in small fixed-size buffers (500 bytes by default), keeping memory overhead flat

## Building

**Requirements:** GCC 15+ (or any C++20-capable compiler), GNU Make

```bash
# Build the executable
make

# Build and run (produces output.wav in the project root)
make run

# Clean build artifacts
make clean
```

The binary is written to `build/audio-editor`.

## How It Works

### Waveforms

Each waveform is a plain function with the signature `double(double timeSeconds, double frequency)`, returning a normalized amplitude in `[-1.0, 1.0]`. Time is always relative to zero rather than the sound's start position — this prevents phase cancellation when combining waveforms.

```
sine     →  sin(2π · frequency · time)
square   →  sign of sine (±1.0)
sawtooth →  linear ramp from -1 to 1, wrapping at each period
triangle →  linear rise and fall, symmetric
```

### Keyframes and Transformations

Amplitude and frequency envelopes are defined as a list of keyframes, each pairing a `value` with a `percent` (position in time from 0.0 to 1.0). Two transformation classes implement the `Frequency` interface:

- **`Linear`** — linearly interpolates between keyframes for amplitude (`g()`), and integrates the area under the curve for frequency (`f()`). The integration step is what makes pitch sweeps sound correct: the wave function receives a cumulative phase value rather than an instantaneous frequency, so a sweep from 440 Hz to 880 Hz produces the expected glide.
- **`Constant`** — returns the same value regardless of position, useful for sounds with no modulation.

### Mixing

The `combine()` function iterates over all active sounds at a given sample index, evaluates each sound's wave function with its modulated frequency and amplitude, sums the results, and divides by the number of active sounds. Sounds that haven't started yet or have already ended are skipped.

### WAV Output

1. A 44-byte header is constructed in memory and written first, encoding sample rate, bit depth, channel count, and data size.
2. The main loop generates samples in `BUFFER_SIZE / bytesPerSample` batches. Each sample is:
   - Combined from all active sounds
   - Quantized from `[-1.0, 1.0]` to a signed 16-bit integer via `quantize<int16_t>()`
   - Converted to little-endian bytes via the `toLittleEndian<T>()` template
   - Written into the current buffer
3. Each full buffer is flushed to the file before the next one begins.

### Template Utilities

```cpp
// Converts a normalized [-1.0, 1.0] sample to a fixed-width integer
template <typename T>
T quantize(double sample, int bitsPerSample);

// Decomposes any integer type into a little-endian byte array
template <typename T>
constexpr std::array<uint8_t, sizeof(T)> toLittleEndian(T value);

// Writes a fixed-size array into a destination buffer at a byte offset
template<typename T, size_t N, size_t M>
void writeAtOffset(const std::array<T, N>& src, size_t offset, std::array<T, M>& dest);
```

## Example

`main.cpp` demonstrates two sine waves with simultaneous frequency and amplitude modulation:

- **Sound 1:** sweeps from 440 Hz to 880 Hz over 2 seconds, with amplitude rising to peak at 50% then falling
- **Sound 2:** sweeps from 880 Hz to 1760 Hz over 2 seconds, with the same amplitude envelope

Both sounds start at the same time and are mixed together. The output is 2 seconds of 44100 Hz, 16-bit mono audio written to `output.wav`.

## Project Structure

```
audio-editor/
├── src/
│   ├── main.cpp      # Entry point — configures and runs the synthesis pipeline
│   ├── wave.cpp      # Wave functions, header creation, mixing, and transformations
│   ├── wave.hpp      # Public interface: Keyframe, Sound, Frequency, Linear, Constant, templates
│   └── buffer.hpp    # Experimental C++20 concepts-based buffer abstraction
├── tests/
│   └── test.cpp      # Placeholder for gtest-based unit tests
├── build/            # Compiled output (generated)
├── makefile
└── README.md
```

## C++20 Features Used

- `std::numbers::pi` for the mathematical constant
- `std::format` for debug output
- `std::span` for non-owning buffer views
- C++20 concepts in `buffer.hpp` (`ByteSink`) for compile-time interface constraints
- Designated initializers for aggregate struct construction

## Reference

WAV file format specification: [soundfile.sapp.org/doc/WaveFormat](https://web.archive.org/web/20250610165025/http://soundfile.sapp.org/doc/WaveFormat/)
