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

void writeHelper(std::span<const std::uint8_t> value, std::size_t offset, std::vector<std::uint8_t>& bytes) {
    if (offset + value.size() > bytes.size()) {
        throw std::out_of_range("writeHelper would write past end of buffer");
    }

    std::copy(value.begin(), value.end(), bytes.begin() + offset);
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

// TODO: should test this more, it could have bugs
// should this function cache it's previous state for performance? too early to tell
// should this function just return the multiplier?


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
        const double normalAmplitude = combine({&s1, &s2}, i);

        // const double x = combine({&s1}, i);
        // const double x = combine({&s2}, i);

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

    std::array<std::uint8_t, 44> header = createHeader(numChannels, sampleRate, bitsPerSample, bytes.size());

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
