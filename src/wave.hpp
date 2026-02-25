#pragma once

#include "buffer.hpp"

#include <vector>
#include <memory>
#include <cassert>
#include <cmath>
#include <array>

struct Keyframe {
    double value;  // should be between 0.0 and 1.0 (but maybe can be negative)
    double percent;     // must be between 0.0 and 1.0
};

// TODO: make all transformations a list so that they can be applied in any order
//          but double check if this is even necessary, we could just need to apply pitch, then amplitude transformatinos

class Frequency {
public:
    // Do I even need this?
    // virtual ~Frequency() = default;
    virtual double f(double percent) = 0;
};


class Linear: public Frequency {
    std::vector<Keyframe> keyframes;
    double cumulativeIntegral(const std::vector<Keyframe> &keyframes, double percent);
    double interpolate(double samplePercent, const std::vector<Keyframe> &keyframes);
public:
    Linear(std::vector<Keyframe> keyframes): keyframes{keyframes} {}
    virtual double f(double percent) {
        return cumulativeIntegral(keyframes, percent);
    };
};

class Constant: public Frequency {
    double value;
public:
    Constant(double value): value{value} {}
    virtual double f(double) {
        return value;
    }
};

struct Sound {
    double (*function)(double, double);
    size_t delaySamples;
    size_t numSamples;
    int sampleRate; // TODO: should this be moved higher up since this can'd differ per function
    // double (*frequency)(double);
    std::unique_ptr<Frequency> frequency;
    std::unique_ptr<Frequency> amplitude;
    // std::vector<Keyframe> linearTransformKeyframes;

    // Sound(double (*function)(size_t, int, double), size_t delaySamples, size_t numSamples, int sampleRate, std::unique_ptr<Frequency> frequency, std::vector<Keyframe> linearTransformKeyframes) 
    //     : function{function}
    //     , delaySamples{delaySamples}
    //     , numSamples{numSamples}
    //     , sampleRate{sampleRate}
    //     , frequency{std::move(frequency)}
    //     , linearTransformKeyframes{linearTransformKeyframes}
    // {}
};


// waves //////////////////////////////////////////////////////////////////////
double squareWave(double timeSeconds, double frequency);
double sawtoothWave(double timeSeconds, double frequency);
double triangleWave(double timeSeconds, double frequency);
double sinWave(double timeSeconds, double frequency);

//
std::array<std::uint8_t, 44> createHeader(std::uint16_t numChannels, std::uint32_t sampleRate, std::uint16_t bitsPerSample, size_t samplesSize);
void createHeader2(IBuffer &buffer, std::uint16_t numChannels, std::uint32_t sampleRate, std::uint16_t bitsPerSample, size_t samplesSize);
double combine(const std::vector<Sound*> &sounds, size_t sampleIndex);

// templates //////////////////////////////////////////////////////////////////
template <typename T>
T quantize(double sample, int bitsPerSample) {
    assert(bitsPerSample % 8 == 0 && bitsPerSample >= 8);
    assert(bitsPerSample / 8 == sizeof(T));

    // TODO: this probably does not need to be computed at runtime
    T maxAmplitude = (std::pow(2, bitsPerSample) / 2) - 1;
    return sample * maxAmplitude;
}

template <typename T>
constexpr std::array<uint8_t, sizeof(T)> toLittleEndian(T value) {
    std::array<uint8_t, sizeof(T)> out{};
    for (size_t i = 0; i < sizeof(T); ++i) {
        out[i] = static_cast<uint8_t>((value >> (8 * i)) & 0xFF);
    }
    return out;
}


