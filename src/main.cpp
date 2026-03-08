#include <iostream>
#include <fstream>
#include <vector>

#include "wave.hpp"
#include "json.hpp"

#define BUFFER_SIZE 500

int main() {
    const int numChannels = 1;
    const int sampleRate = 44100;
    const int bitsPerSample = 16;
    const int bytesPerSample = bitsPerSample / 8;

    const int durationMs = 2000;
    const double frequency = 440.0;
    const size_t numSamples = (durationMs * sampleRate) / 1000;
    const size_t delaySamples = 0;

    std::vector<Keyframe> keyframes = {
        { .value = 0.0, .percent = 0.0 },
        { .value = 1.0, .percent = 0.5 },
        { .value = 0.0, .percent = 1.0 },
    };

    std::vector<Keyframe> amplitudeUpKeyframes = {
        { .value = 0.0, .percent = 0.0 },
        { .value = 1.0, .percent = 1.0 },
    };

    std::vector<Keyframe> amplitudeDownKeyframes = {
        { .value = 1.0, .percent = 0.0 },
        { .value = 0.0, .percent = 1.0 },
    };

    std::vector<Keyframe> s1FrequencyKeyframes = {
        { .value = double(frequency), .percent = 0.0 },
        { .value = double(frequency*2), .percent = 1.0 },
    };

    std::vector<Keyframe> s2FrequencyKeyframes = {
        { .value = double(frequency*2), .percent = 0.0 },
        { .value = double(frequency*4), .percent = 1.0 },
    };

    Sound s1 = {
        .function = sinWave,
        .delaySamples = delaySamples,
        .numSamples = numSamples,
        .sampleRate = sampleRate,
        .frequency = std::make_unique<Linear>(s1FrequencyKeyframes),
        .amplitude = std::make_unique<Linear>(keyframes),
    };

    Sound s2 = {
        .function = sinWave,
        .delaySamples = delaySamples,
        .numSamples = numSamples,
        .sampleRate = sampleRate,
        .frequency = std::make_unique<Linear>(s2FrequencyKeyframes),
        .amplitude = std::make_unique<Linear>(keyframes),
    };

    std::ofstream file("output.wav", std::ios::binary);
    if (!file) {
        std::cerr << "failed to open file\n";
        return 1;
    }

    std::array header = createHeader(numChannels, sampleRate, bitsPerSample, numSamples * bytesPerSample);
    file.write(reinterpret_cast<const char*>(header.data()), header.size());

    for (size_t i = 0; i < numSamples;) {
        std::array<std::uint8_t, BUFFER_SIZE> bytes;

        for (size_t j = 0; j < BUFFER_SIZE / bytesPerSample; j++, i++) {
            const double normalAmplitude = combine({&s1, &s2}, i);
            const auto quantizedAmplitude = quantize<std::int16_t>(normalAmplitude, bitsPerSample);
            const auto byteAmplitude = toLittleEndian(quantizedAmplitude);

            writeAtOffset(byteAmplitude, j * bytesPerSample, bytes);
        }

        file.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    }

    file.close();

    return 0;
}
