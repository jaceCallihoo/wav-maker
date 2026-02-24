#pragma once

#include <vector>
#include <memory>

struct Keyframe {
    double value;  // should be between 0.0 and 1.0 (but maybe can be negative)
    double percent;     // must be between 0.0 and 1.0
};
double keyframeAverage(const std::vector<Keyframe> &keyframes, double percent);

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
public:
    Linear(std::vector<Keyframe> keyframes): keyframes{keyframes} {}
    virtual double f(double percent) {
        return keyframeAverage(keyframes, percent);
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



