#include <arpa/inet.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <array>
#include <span>
#include <cmath>
#include <numbers>
#include <cassert>
#include <algorithm>
#include <format>
#include <algorithm>
#include <memory>

#include "wave.hpp"

#define assertm(exp, msg) assert((void(msg), exp))

// TODO: shoulde the function names be camel case? Check google formatting standards
// TODO: check the class field names as well


template <typename T>
constexpr std::array<uint8_t, sizeof(T)> toLittleEndian(T value) {
    std::array<uint8_t, sizeof(T)> out{};
    for (size_t i = 0; i < sizeof(T); ++i) {
        out[i] = static_cast<uint8_t>((value >> (8 * i)) & 0xFF);
    }
    return out;
}

void writeHelper(std::span<const std::uint8_t> value, std::size_t offset, std::vector<std::uint8_t>& bytes) {
    if (offset + value.size() > bytes.size()) {
        throw std::out_of_range("writeHelper would write past end of buffer");
    }

    std::copy(value.begin(), value.end(), bytes.begin() + offset);
}

std::vector<std::uint8_t> createHeader(std::uint16_t numChannels, std::uint32_t sampleRate, std::uint16_t bitsPerSample, size_t samplesSize) {
    const static std::size_t offsetChunkId = 0;
    const static std::size_t offsetChunkSize = 4;
    const static std::size_t offsetFormat = 8;
    const static std::size_t offsetSubchunk1Id = 12;
    const static std::size_t offsetSubchunk1Size = 16;
    const static std::size_t offsetAudioFormat = 20;
    const static std::size_t offsetNumChannels = 22;
    const static std::size_t offsetSampleRate = 24;
    const static std::size_t offsetByteRate = 28;
    const static std::size_t offsetBlockAlign = 32;
    const static std::size_t offsetBitsPerSample = 34;
    const static std::size_t offsetSubchunk2Id = 36;
    const static std::size_t offsetSubchunk2Size = 40;
    // const static std::size_t offsetData = 44;

    const static constexpr std::array<std::uint8_t, 4> chunkId = {'R', 'I', 'F', 'F'};
    const static constexpr std::array<std::uint8_t, 4> format = {'W', 'A', 'V', 'E'};
    const static constexpr std::array<std::uint8_t, 4> subchunk1Id = {'f', 'm', 't', ' '};
    const static constexpr std::array<std::uint8_t, 4> subchunk1Size = toLittleEndian<std::uint32_t>(16);
    const static constexpr std::array<std::uint8_t, 2> audioFormat = toLittleEndian<std::uint16_t>(1);
    const static constexpr std::array<std::uint8_t, 4> subchunk2Id = {'d', 'a', 't', 'a'};

    std::array<std::uint8_t, 4> chunkSize = toLittleEndian<std::uint32_t>(36 + samplesSize);
    std::array<std::uint8_t, 4> byteRate = toLittleEndian<std::uint32_t>(sampleRate * numChannels * bitsPerSample / 8);
    std::array<std::uint8_t, 2> blockAlign = toLittleEndian<std::uint16_t>(numChannels * bitsPerSample / 8);
    std::array<std::uint8_t, 4> subchunk2Size = toLittleEndian<std::uint32_t>(samplesSize);

    std::vector<std::uint8_t> bytes(44);

    // TODO: should this take a writer as an argument insteado of writing to a buffer?
    writeHelper(chunkId, offsetChunkId, bytes);
    writeHelper(chunkSize, offsetChunkSize, bytes);
    writeHelper(format, offsetFormat, bytes);
    writeHelper(subchunk1Id, offsetSubchunk1Id, bytes);
    writeHelper(subchunk1Size, offsetSubchunk1Size, bytes);
    writeHelper(audioFormat, offsetAudioFormat, bytes);
    writeHelper(toLittleEndian<std::uint16_t>(numChannels), offsetNumChannels, bytes);
    writeHelper(toLittleEndian<std::uint32_t>(sampleRate), offsetSampleRate, bytes);
    writeHelper(byteRate, offsetByteRate, bytes);
    writeHelper(blockAlign, offsetBlockAlign, bytes);
    writeHelper(toLittleEndian<std::uint16_t>(bitsPerSample), offsetBitsPerSample, bytes);
    writeHelper(subchunk2Id, offsetSubchunk2Id, bytes);
    writeHelper(subchunk2Size, offsetSubchunk2Size, bytes);

    return bytes;
}

// bitsPerSample: how many bits are used to record a sample (not important until later when writing)
// apmlitude: distance of the wave sample from 0 (also not super important to the calc)

// sampleRate: how many samples of the wave are there in a seccond
// frequency: number of sound wave cycels per second
//      1Hz = 1 cycle per second
// wavelength: 
// waves per second? -> how many waves have been before this point?
//      frequency 

// i: the current sample
// double(i) / sampleRate: position within the second ... is this even useful? don't we want position within the wave?
//      maybe this is useful when combined with frequency
// period: 1 / frequency

// TODO: can this returen an array somehow instead of a vector?
// TODO: make a generator version of this function
// TODO: make a stateless version of this function that gets a saple at a given index

// TODO: maybe make this later, or find a way to modify the square wave maybe? Maybe have a translation
// template <typename T>
// std::vector<T> pulseWave(int durationMs, int sampleRate, int bitsPerSample, int frequency, T amplitude) {
// }

double sawtoothWave(double timeSeconds, double frequency) {
    assert(timeSeconds >= 0);
    assert(frequency > 0);

    const double periodsElapsed = timeSeconds * frequency;
    const double partialPeriod = periodsElapsed - floor(periodsElapsed);

    if (partialPeriod < 0.5) { 
        return 2 * partialPeriod;
    } else {
        return (2 * partialPeriod) + -2;
    }
}

double triangleWave(double timeSeconds, double frequency) {
    assert(timeSeconds >= 0);
    assert(frequency > 0);

    const double periodsElapsed = timeSeconds * frequency;
    const double partialPeriod = periodsElapsed - floor(periodsElapsed);

    if (partialPeriod < 0.25) { 
        return 4 * partialPeriod;
    } else if (partialPeriod < 0.75) {
        return (-4 * partialPeriod) + 2;
    } else {
        return (4 * partialPeriod) + -4;
    }
}

double sinWave(double timeSeconds, double frequency) {
    assert(timeSeconds >= 0);
    assert(frequency > 0);

    return std::sin(2 * std::numbers::pi * frequency * timeSeconds);
}

// assumes sorted data
double keyframeAverage(std::vector<Keyframe> keyframes) {
    assert(keyframes.size() >= 2);
    assert(keyframes[0].percent == 0.0);
    assert(keyframes[keyframes.size()-1].percent == 1.0);

    double result = 0.0;

    for (size_t i = 1; i < keyframes.size(); i++) {
        const double gapSize = (keyframes[i].percent - keyframes[i-1].percent);
        const double averageMultiplier = (keyframes[i].value + keyframes[i-1].value) / 2;

        result += gapSize * averageMultiplier;
    }

    return result;
}

// TODO: should test this more, it could have bugs
double keyframeAverage(const std::vector<Keyframe> &keyframes, double percent) {
    assert(percent >= 0);
    assert(percent <= 1.0);
    assert(keyframes.size() >= 2);
    assert(keyframes[0].percent == 0.0);
    assert(keyframes[keyframes.size()-1].percent == 1.0);

    if (percent == 0.0) {
        return keyframes[0].value;
    }


    double result = 0.0;

    // add up all the previous gaps
    size_t i = 1;
    while (i < keyframes.size() && keyframes[i].percent < percent) {
        const double averageMultiplier = (keyframes[i].value + keyframes[i-1].value) / 2;
        const double gapSize = (keyframes[i].percent - keyframes[i-1].percent);

        result += gapSize * averageMultiplier;

        i++;
    }

    if (i < keyframes.size()) {
        double keyframeMultiplierDifference = keyframes[i-1].value - keyframes[i].value;  // rise
        double keyframePercentDifference = keyframes[i-1].percent - keyframes[i].percent;           // run
        double slope = keyframeMultiplierDifference / keyframePercentDifference;
        double x = percent - keyframes[i-1].percent;
        double y = x * slope + keyframes[i-1].value;

        result += ((y + keyframes[i-1].value) / 2) * x;

    }

    result /= percent;

    return result;
}

// double squareWave2(size_t sampleIndex, int sampleRate, const std::vector<Keyframe> &amplitudeUpKeyframes, double percent) {
//     assert(sampleRate > 0);
//     assert(amplitudeUpKeyframes.size() > 0);
// 
//     X x = keyframeAverage(amplitudeUpKeyframes, percent);
// 
//     const double timeSeconds = double(sampleIndex) / sampleRate;
//     const double periodsElapsed = timeSeconds * x.previousAverage;
//     const double partialPeriod = periodsElapsed - floor(periodsElapsed);
// 
//     if (partialPeriod < 0.5) {
//         return 1.0;
//     } else {
//         return -1.0;
//     }
// }

double squareWave(double timeSeconds, double frequency) {
    assert(timeSeconds >= 0);
    assert(frequency > 0);

    const double periodsElapsed = timeSeconds * frequency;
    const double partialPeriod = periodsElapsed - floor(periodsElapsed);

    if (partialPeriod < 0.5) {
        return 1.0;
    } else {
        return -1.0;
    }
}

template <typename T>
T quantize(double sample, int bitsPerSample) {
    assert(bitsPerSample % 8 == 0 && bitsPerSample >= 8);
    assert(bitsPerSample / 8 == sizeof(T));

    // TODO: this probably does not need to be computed at runtime
    T maxAmplitude = (std::pow(2, bitsPerSample) / 2) - 1;
    return sample * maxAmplitude;
}

// vertical transforms ////////////////////////////////////////////////////////
// TODO: probably don't need these anymore
void transformScale(std::vector<double> &samples, double scaler) {
    for (auto &sample : samples) {
        sample *= scaler;
    }
}

// TODO: if this is still used it can just be a scale with -1.0
void transformInvert(std::vector<double> &samples) {
    for (auto &sample : samples) {
        sample = -sample;
    }
}

void transformOffset(std::vector<double> &samples, double offset) {
    for (auto &sample : samples) {
        sample += offset;
    }
}

// void transformBezier(std::vector<double> &samples, double x1, double y1, double x2, double y2) {
//     // TODO: maybe normalize first?
//     for (auto &sample : samples) {
//         sample = 0;
//     }
// }

// should this function cache it's previous state for performance? too early to tell
// should this function just return the multiplier?
double transformLinear(double samplePercent, const std::vector<Keyframe> &keyframes) {
    assert(samplePercent >= 0);
    assert(samplePercent <= 1.0);
    assert(keyframes.size() > 0);

    // hold first / hold last

    if (keyframes.size() == 1) {
        return keyframes[0].value;
    }

    // if we're before the first keyframe
    // return sample * keframe
    if (samplePercent <= keyframes[0].percent) {
        return keyframes[0].value;
    }

    // if we're after the last keyframe
    if (samplePercent >= keyframes[keyframes.size() - 1].percent) {
        return keyframes[keyframes.size() - 1].value;
    }

    size_t leftKeyframeIndex = 0;
    size_t rightKeyframeIndex = 1;
    while (samplePercent >= keyframes[rightKeyframeIndex].percent) {
        leftKeyframeIndex++;
        rightKeyframeIndex++;
    }

    double keyframePercentDifference = keyframes[leftKeyframeIndex].percent - keyframes[rightKeyframeIndex].percent;
    double keyframeMultiplierDifference = keyframes[leftKeyframeIndex].value - keyframes[rightKeyframeIndex].value;
    double slope = keyframeMultiplierDifference / keyframePercentDifference;
    double sampleMultiplier = (samplePercent - keyframes[leftKeyframeIndex].percent) * slope + keyframes[leftKeyframeIndex].value;

    return sampleMultiplier;
}

double combine2(const std::vector<Sound*> &sounds, size_t sampleIndex) {
    double combined = 0;

    for (const auto &sound : sounds) {
        const auto& function = sound->function;
        const auto& delaySamples = sound->delaySamples;
        const auto& numSamples = sound->numSamples;
        const auto& sampleRate = sound->sampleRate;
        const auto& frequency = sound->frequency;
        const auto& amplitude = sound->amplitude;

        // const auto& linearTransformKeyframes = sound->linearTransformKeyframes;

        if (sampleIndex < delaySamples) {
            continue;
        }

        if (sampleIndex >= delaySamples + numSamples) {
            continue;
        }

        const double timeSeconds = double(sampleIndex) / sampleRate;
        const double percent = double(sampleIndex - delaySamples) / numSamples;

        double current = function(timeSeconds, frequency->f(percent));

        // const double activeSamplePercent = double(sampleIndex - delaySamples) / numSamples;

        // if (linearTransformKeyframes.size() > 0) {
        //     combined *= transformLinear(activeSamplePercent, linearTransformKeyframes);
        // }

        current *= amplitude->f(percent);

        combined += current;
    }

    combined /= sounds.size();

    return combined;
}

int main() {
    
    // set input values
    const int numChannels = 1;
    const int sampleRate = 44100;
    const int bitsPerSample = 16;

    // create sound with input values
    // std::vector<std::uint8_t> data = createSound(numChannels, sampleRate, bitsPerSample);
    // create wav file with input values
    const int durationMs = 2000;
    const double frequency = 440.0;
    // const std::int16_t amplitude = 16383;
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

    // const double frequency = 440.0;

    std::vector<Keyframe> s1FrequencyKeyframes = {
        { .value = double(frequency), .percent = 0.0 },
        { .value = double(frequency*2), .percent = 1.0 },
    };

    std::vector<Keyframe> s2FrequencyKeyframes = {
        { .value = double(frequency*2), .percent = 0.0 },
        { .value = double(frequency*4), .percent = 1.0 },
    };

    // for (double percent = 0.0; percent <= 1.0; percent += 0.01) {
    //     X x = keyframeAverage(amplitudeUpKeyframes, percent);
    //     std::cout << std::format("percent: {}\npreviousAverage: {}\ncurrent: {}\n\n", percent, x.previousAverage, x.current);
    // }
    // exit(0);

    // Sound s1(sinWave, delaySamples, numSamples, sampleRate, std::make_unique<Constant>(440), keyframes);
    // Sound s1 = {
    //     .function = sinWave,
    //     .delaySamples = 0,
    //     .numSamples = numSamples,
    //     .sampleRate = sampleRate,
    //     .frequency = std::make_unique<Constant>(frequency),
    //     .linearTransformKeyframes = keyframes,
    // };

    Sound s1 = {
        .function = sinWave,
        .delaySamples = 0,
        .numSamples = numSamples,
        .sampleRate = sampleRate,
        .frequency = std::make_unique<Linear>(s1FrequencyKeyframes),
        .amplitude = std::make_unique<Linear>(keyframes),
    };

    Sound s2 = {
        .function = sinWave,
        .delaySamples = 0,
        .numSamples = numSamples,
        .sampleRate = sampleRate,
        .frequency = std::make_unique<Linear>(s2FrequencyKeyframes),
        .amplitude = std::make_unique<Linear>(keyframes),
    };


    // TEST //////////////
    // const size_t i = 699;
    // const double percent = double(i) / numSamples;
    // X x = keyframeAverage(amplitudeUpKeyframes, percent);
    // std::cout << std::format("outside x: current = {}, previousAverage = {}\n", x.current, x.previousAverage);
    // const double normalAmplitude2 = squareWave(i, sampleRate, x.previousAverage);
    // const double normalAmplitude = squareWave2(i, sampleRate, amplitudeUpKeyframes, double(i)/ numSamples);
    // std::cout << std::format("diff: {}\na: {}\nb: {}\npercent: {}\n\n", i, normalAmplitude2, normalAmplitude, double(i)/numSamples);
    // exit(0);

    ///////////////



    // const int UP = 1;
    // const int DOWN = 2;

    // int state = UP;
    // size_t previousPeriodStart = 0;

    std::vector<std::uint8_t> bytes(numSamples * (bitsPerSample / 8));
    for (size_t i = 0; i < numSamples; i++) {
        const double normalAmplitude = combine2({&s1, &s2}, i);

        // const double x = combine2({&s1}, i);
        // const double x = combine2({&s2}, i);

        // if (x > 0.0 && state == DOWN) {
        //     state = UP;
        //     size_t numSamplesInPeriod = i - previousPeriodStart;
        //     double periodDuration = double(numSamplesInPeriod) / sampleRate;
        //     double frequency = 1 / periodDuration;
        //     std::cout << std::format("wave frequency: {}\nperiodDuration: {}s\nnumSamplesInPeriod: {}\n\n", frequency, periodDuration, numSamplesInPeriod);
        //     previousPeriodStart = i;
        // } else if (x < 0.0 && state == UP) {
        //     state = DOWN;
        //     // do nothing else here since it's only half
        // }

        // std::string buffer(20, ' ');
        // size_t index = std::min((1 + normalAmplitude) * 10, double(19));
        // buffer[index] = '.';
        // std::cout << buffer << '\n';

        const auto quantizedAmplitude = quantize<std::int16_t>(normalAmplitude, bitsPerSample);
        const auto byteAmplitude = toLittleEndian(quantizedAmplitude);
        writeHelper(byteAmplitude, i * (bitsPerSample / 8), bytes);
        // bytes.insert(bytes.end(), byteAmplitude.begin(), byteAmplitude.end());
    }

    // for (size_t i = 0; i < bytes.size(); i++) {
    //     std::cout << std::hex << int(bytes[i]) << " ";
    // }
    // std::cout << std::endl;

    // Wav wave(numChannels, sampleRate, bitsPerSample, bytes);

    // std::vector<std::uint8_t> buffer = wave.toBytes();

    std::vector<std::uint8_t> header = createHeader(numChannels, sampleRate, bitsPerSample, bytes.size());

    //////////////////////////////////////////////////////////////////////
    // write to file

    std::ofstream file("output.wav", std::ios::binary);
    if (!file) {
        std::cerr << "failed to open file\n";
        return 1;
    }


    // file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    file.write(reinterpret_cast<const char*>(header.data()), header.size());
    file.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    // optional, this is called when the file leaves scope
    file.close();

    return 0;
}
