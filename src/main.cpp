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

// TODO: fix this so the input does not need to be converted into a span
template <typename T>
std::vector<std::uint8_t> allToLittleEndian(std::span<T> values) {
    std::vector<std::uint8_t> bytes(values.size() * sizeof(T));

    for (std::size_t i = 0; i < values.size(); i++) {
        size_t offset = i * sizeof(T);
        writeHelper(toLittleEndian(values[i]), offset, bytes);
    }

    return bytes;
}


class Wav {
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
    const static std::size_t offsetData = 44;

    const static constexpr std::array<std::uint8_t, 4> chunkId = {'R', 'I', 'F', 'F'};
    const static constexpr std::array<std::uint8_t, 4> format = {'W', 'A', 'V', 'E'};
    const static constexpr std::array<std::uint8_t, 4> subchunk1Id = {'f', 'm', 't', ' '};
    const static constexpr std::array<std::uint8_t, 4> subchunk1Size = toLittleEndian<std::uint32_t>(16);
    const static constexpr std::array<std::uint8_t, 2> audioFormat = toLittleEndian<std::uint16_t>(1);
    const static constexpr std::array<std::uint8_t, 4> subchunk2Id = {'d', 'a', 't', 'a'};

    std::uint16_t numChannels;
    std::uint32_t sampleRate;
    std::uint16_t bitsPerSample;
    // TODO: maybe this should be a pointer? Or have some sort of memory optimization for large files
    std::vector<std::uint8_t> data;

    std::uint32_t deriveChunkSize() { return 36 + deriveSubchunk2Size(); }
    std::uint32_t deriveByteRate() { return sampleRate * numChannels * bitsPerSample / 8; }
    std::uint16_t deriveBlockAlign() { return numChannels * bitsPerSample / 8; }
    std::uint32_t deriveSubchunk2Size() { return data.size(); }

    size_t getTotalSize() { return deriveChunkSize() + 8; }

public:
    Wav(std::uint16_t numChannels, std::uint32_t sampleRate, std::uint16_t bitsPerSample, std::vector<std::uint8_t> data)
        : numChannels(numChannels), sampleRate(sampleRate), bitsPerSample(bitsPerSample), data(data)
    {}

    std::vector<std::uint8_t> toBytes() {
        std::vector<std::uint8_t> bytes(getTotalSize());

        std::array<std::uint8_t, 2> serializedNumChannels = toLittleEndian(numChannels);
        std::array<std::uint8_t, 4> serializedSampleRate = toLittleEndian(sampleRate);
        std::array<std::uint8_t, 2> serializedBitsPerSample = toLittleEndian(bitsPerSample);
        std::array<std::uint8_t, 4> serializedChunkSize = toLittleEndian(deriveChunkSize());
        std::array<std::uint8_t, 4> serializedByteRate = toLittleEndian(deriveByteRate());
        std::array<std::uint8_t, 2> serializedBlockAlign = toLittleEndian(deriveBlockAlign());
        std::array<std::uint8_t, 4> serializedSubchunk2Size = toLittleEndian(deriveSubchunk2Size());

        writeHelper(chunkId, offsetChunkId, bytes);
        writeHelper(serializedChunkSize, offsetChunkSize, bytes);
        writeHelper(format, offsetFormat, bytes);
        writeHelper(subchunk1Id, offsetSubchunk1Id, bytes);
        writeHelper(subchunk1Size, offsetSubchunk1Size, bytes);
        writeHelper(audioFormat, offsetAudioFormat, bytes);
        writeHelper(serializedNumChannels, offsetNumChannels, bytes);
        writeHelper(serializedSampleRate, offsetSampleRate, bytes);
        writeHelper(serializedByteRate, offsetByteRate, bytes);
        writeHelper(serializedBlockAlign, offsetBlockAlign, bytes);
        writeHelper(serializedBitsPerSample, offsetBitsPerSample, bytes);
        writeHelper(subchunk2Id, offsetSubchunk2Id, bytes);
        writeHelper(serializedSubchunk2Size, offsetSubchunk2Size, bytes);
        writeHelper(data, offsetData, bytes);

        return bytes;
    }
};

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

double sawtoothWave2(size_t sampleIndex, int sampleRate, int frequency) {
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double timeSeconds = double(sampleIndex) / sampleRate;
    const double periodsElapsed = timeSeconds * frequency;
    const double partialPeriod = periodsElapsed - floor(periodsElapsed);

    if (partialPeriod < 0.5) { 
        return 2 * partialPeriod;
    } else {
        return (2 * partialPeriod) + -2;
    }
}

std::vector<double> sawtoothWave(size_t numSamples, int sampleRate, int frequency) {
    assert(numSamples > 0);
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double period = double(1) / frequency;

    std::vector<double> samples(numSamples);

    for (std::size_t i = 0; i < numSamples; i++) {
        const double timeSeconds = double(i) / sampleRate;
        const double periodTimeSeconds = std::fmod(timeSeconds, period);
        const double periodPercent = periodTimeSeconds / period;
        if (periodPercent < 0.5) { 
            samples[i] = 2 * periodPercent;
        } else {
            samples[i] = (2 * periodPercent) + -2;
        }
    }

    return samples;
}

double triangleWave2(size_t sampleIndex, int sampleRate, int frequency) {
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double timeSeconds = double(sampleIndex) / sampleRate;
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

std::vector<double> triangleWave(size_t numSamples, int sampleRate, int frequency) {
    assert(numSamples > 0);
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double period = double(1) / frequency;

    std::vector<double> samples(numSamples);

    for (std::size_t i = 0; i < numSamples; i++) {
        const double timeSeconds = double(i) / sampleRate;
        const double periodTimeSeconds = std::fmod(timeSeconds, period);
        const double periodPercent = periodTimeSeconds / period;
        if (periodPercent < 0.25) { 
            samples[i] = 4 * periodPercent;
        } else if (periodPercent < 0.75) {
            samples[i] = (-4 * periodPercent) + 2;
        } else {
            samples[i] = (4 * periodPercent) + -4;
        }
    }

    return samples;
}

double sinWave2(size_t sampleIndex, int sampleRate, int frequency) {
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double timeSeconds = double(sampleIndex) / sampleRate;

    return std::sin(2 * std::numbers::pi * frequency * timeSeconds);
}

std::vector<double> sinWave(size_t numSamples, int sampleRate, int frequency) {
    assert(numSamples > 0);
    assert(sampleRate > 0);
    assert(frequency > 0);

    std::vector<double> samples(numSamples);

    for (std::size_t i = 0; i < numSamples; i++) {
        const double timeSeconds = double(i) / sampleRate;
        samples[i] = std::sin(2 * std::numbers::pi * frequency * timeSeconds);
    }

    return samples;
}

double squareWave2(size_t sampleIndex, int sampleRate, int frequency) {
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double timeSeconds = double(sampleIndex) / sampleRate;
    const double periodsElapsed = timeSeconds * frequency;
    const double partialPeriod = periodsElapsed - floor(periodsElapsed);    // this should better distribute the error

    if (partialPeriod < 0.5) {
        return 1.0;
    } else {
        return -1.0;
    }
}

std::vector<double> squareWave(size_t numSamples, int sampleRate, int frequency) {
    assert(numSamples > 0);
    assert(sampleRate > 0);
    assert(frequency > 0);

    const double period = double(1) / frequency;

    std::vector<double> samples(numSamples);

    for (std::size_t i = 0; i < numSamples; i++) {
        const double timeSeconds = double(i) / sampleRate;
        const double periodTimeSeconds = std::fmod(timeSeconds, period);
        const double periodPercent = periodTimeSeconds / period;
        if (periodPercent < 0.5) {
            samples[i] = 1.0;
        } else {
            samples[i] = -1.0;
        }
    }

    return samples;
}

template <typename T>
T quantize2(double sample, int bitsPerSample) {
    assert(bitsPerSample % 8 == 0 && bitsPerSample >= 8);
    assert(bitsPerSample / 8 == sizeof(T));

    // TODO: this probably does not need to be computed at runtime
    T maxAmplitude = (std::pow(2, bitsPerSample) / 2) - 1;
    return sample * maxAmplitude;
}


// TODO: see if you can generalize the vector
// TODO: maybe ranges wolud be better than boomer loop
template <typename T>
std::vector<T> quantize(std::vector<double> samples, int bitsPerSample) {
    assert(bitsPerSample % 8 == 0 && bitsPerSample >= 8);
    assert(bitsPerSample / 8 == sizeof(T));

    T maxAmplitude = (std::pow(2, bitsPerSample) / 2) - 1;
    std::vector<T> result(samples.size());
    for (size_t i = 0; i < samples.size(); i++) {
        result[i] = samples[i] * maxAmplitude;
    }

    return result;
}


// vertical transforms ////////////////////////////////////////////////////////
void transformScale(std::vector<double> &samples, double scaler) {
    for (auto &sample : samples) {
        sample *= scaler;
    }
}

void transformInvert(std::vector<double> &samples) {
    for (auto &sample : samples) {
        sample = -sample;
    }
}

void transformInvert(std::vector<double> &samples, double offset) {
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

struct Keyframe {
    double multiplier;  // should be between 0.0 and 1.0 (but maybe can be negative)
    double percent;     // must be between 0.0 and 1.0
};

// TODO: add start and end percent
// assumes that keyframes are in order
void transformLinear(std::vector<double> &samples, const std::vector<Keyframe> &keyframes) {
    // do I need to add 0% and 1.00% (if not there) to make the algo work?
    //      and if so, what should they be (0, 1.0, based of of the next/previous 2 keyframes)

    // TODO: maybe this is an if check that just returns
    // or maybe keyframe 0% and keyframe 100% should be added before assigning leftKeyframe and rightKeyframe
    assert(keyframes.size() > 0);

    // TODO: maybe convert the keyframe percents into indicies? (maybe this is already happening idk)

    // auto nextKeyframePercent = keyframes.begin();
    size_t numSamples = samples.size();
    double currentPercent = 0.0;
    size_t currentSampleIndex = 0;
    size_t leftKeyframeIndex = 0;
    size_t rightKeyframeIndex = 1;
    for (auto &sample : samples) {
        // find surrounding keyframes
        // TODO: can this seg fault?
        while (currentPercent >= keyframes[rightKeyframeIndex].percent) {
            leftKeyframeIndex++;
            rightKeyframeIndex++;
        }

        // interpolate
        double keyframePercentDifference = keyframes[leftKeyframeIndex].percent - keyframes[rightKeyframeIndex].percent;
        double keyframeMultiplierDifference = keyframes[leftKeyframeIndex].multiplier - keyframes[rightKeyframeIndex].multiplier;
        double slope = keyframeMultiplierDifference / keyframePercentDifference;
        double sampleMultiplier = (currentPercent - keyframes[leftKeyframeIndex].percent) * slope + keyframes[leftKeyframeIndex].multiplier;
        sample *= sampleMultiplier;

        // update sample position
        currentSampleIndex++;
        currentPercent = double(currentSampleIndex) / numSamples;
    }
}

///////////////////////////////////////////////////////////////////////////////


// template<typename Vectors>
// void vizualize(int bufferSize, const Vectors&... vectors) {
//     assert(bufferSize > 0);
//     static_assert((std::is_same_v<Vectors, std::vector<double>> && ...),
//                   "All arguments after int must be std::vector<double>");
// 
//     for (const auto &sample : samples) {
//         std::string buffer(20, ' ');
//         size_t index = std::min((1 + sample) * 10, double(19));
//         buffer[index] = '.';
//         std::cout << buffer << '\n';
//     }
// }

struct Sound {
    int sampleRate;
    size_t delaySamples;
    std::vector<double> samples;
};

struct Generator {
    double (*function)(size_t, int, int);
    size_t delaySamples;
    size_t numSamples;
    int sampleRate; // TODO: should this be moved higher up since this can'd differ per function
    int frequency;
};

// should this return a double?
double combine2(std::vector<Generator> generators, size_t sampleIndex) {
    double combined = 0;

    for (const auto &generator : generators) {
        const auto& function = generator.function;
        const auto& delaySamples = generator.delaySamples;
        const auto& numSamples = generator.numSamples;
        const auto& sampleRate = generator.sampleRate;
        const auto& frequency = generator.frequency;

        if (sampleIndex < delaySamples) {
            continue;
        }

        if (sampleIndex > delaySamples + numSamples) {
            continue;
        }

        combined += function(sampleIndex, sampleRate, frequency) / generators.size();
    }

    combined /= generators.size();

    return combined;
}

Sound combine(const std::vector<Sound> &sounds) {
    assert(sounds.size() > 0);
    for (size_t i = 1; i < sounds.size(); i++) {
        assert(sounds[i].sampleRate == sounds[i-1].sampleRate);
    }
    // assert all sampleRate are the same 
    // can we do this at compile time by making sampleRate a type?

    auto getLastSampleIndex = [](Sound sound){
        return sound.delaySamples + sound.samples.size() - 1;
    };

    size_t earliestSampleIndex = sounds[0].delaySamples;
    size_t latestSampleIndex = getLastSampleIndex(sounds[0]);
    for (size_t i = 1; i < sounds.size(); i++) {
        earliestSampleIndex = std::min(earliestSampleIndex, sounds[i].delaySamples);
        latestSampleIndex = std::max(latestSampleIndex, getLastSampleIndex(sounds[i]));
    }

    const size_t resultSize = latestSampleIndex - earliestSampleIndex + 1;
    Sound result = {
        .sampleRate = sounds[0].sampleRate,
        .delaySamples = earliestSampleIndex,
        .samples = std::vector<double>(resultSize),
    };

    for (size_t i = 0; i < sounds.size(); i++) {
        const size_t relativeSoundStartIndex = sounds[i].delaySamples - result.delaySamples;
        for (size_t j = 0; j < sounds[i].samples.size(); j++) {
            result.samples[relativeSoundStartIndex + j] += sounds[i].samples[j];
        }
    }

    for (size_t i = 0; i < result.samples.size(); i++) {
        result.samples[i] /= sounds.size();
    }

    return result;
}

void crop(Sound &sound, size_t leftCrop, size_t rightCrop) {

    const size_t delaySamples = sound.delaySamples;
    const auto& samples = sound.samples;

    // adjust the delay
    sound.delaySamples = std::max<size_t>(0, delaySamples - leftCrop);

    sound.samples.erase(samples.begin(), samples.begin() + std::clamp<size_t>(leftCrop - delaySamples, 0, samples.size()));
    // Idk if I need a bounds check here
    sound.samples.resize(samples.size() - rightCrop);
}

int main() {
    std::cout << "hello world" << std::endl;

    //////////////////////////////////////////////////////////////////////
    
    // set input values
    const int numChannels = 1;
    const int sampleRate = 44100;
    const int bitsPerSample = 16;

    // create sound with input values
    // std::vector<std::uint8_t> data = createSound(numChannels, sampleRate, bitsPerSample);
    // create wav file with input values
    const int durationMs = 2000;
    const int frequency = 440;
    // const std::int16_t amplitude = 16383;
    const size_t numSamples = (durationMs * sampleRate) / 1000;

    Generator g1 = {
        .function = sinWave2,
        .delaySamples = 50000,
        .numSamples = numSamples,
        .sampleRate = sampleRate,
        .frequency = frequency,
    };

    std::vector<std::uint8_t> bytes;
    for (size_t i = 0; i < numSamples; i++) {
        const double normalAmplitude = combine2(std::vector<Generator>({g1}), i);

        std::string buffer(20, ' ');
        size_t index = std::min((1 + normalAmplitude) * 10, double(19));
        buffer[index] = '.';
        std::cout << buffer << '\n';

        const auto quantizedAmplitude = quantize2<std::int16_t>(normalAmplitude, bitsPerSample);
        const auto byteAmplitude = toLittleEndian(quantizedAmplitude);
        bytes.insert(bytes.end(), byteAmplitude.begin(), byteAmplitude.end());
    }

    // std::vector<std::int16_t> squareSamples = squareWave(durationMs, sampleRate, bitsPerSample, frequency, amplitude);
    // std::vector<std::int16_t> sinSamples = sinWave(durationMs, sampleRate, bitsPerSample, frequency, amplitude);

    // std::vector<std::int16_t> samples(squareSamples.size());

    // std::transform(squareSamples.begin(), squareSamples.end(), sinSamples.begin(), samples.begin(), [](std::int16_t x, std::int16_t y) {
    //             return (x / 2 + y / 2) + ((x % 2 + y % 2) / 2);
    //         });

    std::vector<double> samples2 = sinWave(numSamples, sampleRate, frequency);
    // for (int i = 0; i < 100; i++) {
    //     std::cout << std::format("{}\n", samples2[i]);
    // }

    std::vector<Keyframe> keyframes = {
        { .multiplier = 0.0, .percent = 0.0},
        { .multiplier = 1.0, .percent = 0.5},
        { .multiplier = 0.0, .percent = 1.0},
    };
    // transformLinear(samples2, keyframes);
    // vizualize(samples2);
    Sound s = {
        .sampleRate = sampleRate,
        .delaySamples = 0,
        .samples = samples2,
    };

    std::vector<double> samples3(samples2);
    transformInvert(samples3);
    Sound s3 = {
        .sampleRate = sampleRate,
        .delaySamples = 188200,
        .samples = samples3,
    };

    Sound s2 = combine(std::vector<Sound>({s, s3}));
    // crop(s2, 88200, 0);
    std::cout << std::format("samples2.size() {} - s2.samples.size() {}\n", samples2.size(), s2.samples.size());

    std::vector<std::int16_t> samples = quantize<std::int16_t>(s2.samples, bitsPerSample);

    // for (int i = 0; i < 100; i++) {
    //     std::cout << std::format("{}\n", samples[i]);
    // }
    // std::cout << std::endl;

    std::vector<std::uint8_t> data = allToLittleEndian(std::span<std::int16_t>{samples});

    Wav wave(numChannels, sampleRate, bitsPerSample, bytes);

    std::vector<std::uint8_t> buffer = wave.toBytes();

    //////////////////////////////////////////////////////////////////////
    // write to file

    std::ofstream file("output.wav", std::ios::binary);
    if (!file) {
        std::cerr << "failed to open file\n";
        return 1;
    }

    file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    // optional, this is called when the file leaves scope
    file.close();

    return 0;
}
