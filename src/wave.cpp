#include "wave.hpp"

#include <cassert>
#include <cmath>
// TODO: remove
#include <span>

// Linear /////////////////////////////////////////////////////////////////////

double Linear::cumulativeIntegral(const std::vector<Keyframe> &keyframes, double percent) {
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

double Linear::interpolate(double samplePercent, const std::vector<Keyframe> &keyframes) {
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

// waves //////////////////////////////////////////////////////////////////////
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


// 
void createHeader2(IBuffer &buffer, std::uint16_t numChannels, std::uint32_t sampleRate, std::uint16_t bitsPerSample, size_t samplesSize) {
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

    std::array<std::uint8_t, 44> bytes;


    // TODO: should this take a writer as an argument insteado of writing to a buffer?
    buffer.write(chunkId, chunkId.size());
    // writeHelper(chunkId, offsetChunkId, bytes);
    // writeHelper(chunkSize, offsetChunkSize, bytes);
    // writeHelper(format, offsetFormat, bytes);
    // writeHelper(subchunk1Id, offsetSubchunk1Id, bytes);
    // writeHelper(subchunk1Size, offsetSubchunk1Size, bytes);
    // writeHelper(audioFormat, offsetAudioFormat, bytes);
    // writeHelper(toLittleEndian<std::uint16_t>(numChannels), offsetNumChannels, bytes);
    // writeHelper(toLittleEndian<std::uint32_t>(sampleRate), offsetSampleRate, bytes);
    // writeHelper(byteRate, offsetByteRate, bytes);
    // writeHelper(blockAlign, offsetBlockAlign, bytes);
    // writeHelper(toLittleEndian<std::uint16_t>(bitsPerSample), offsetBitsPerSample, bytes);
    // writeHelper(subchunk2Id, offsetSubchunk2Id, bytes);
    // writeHelper(subchunk2Size, offsetSubchunk2Size, bytes);

    // return bytes;
}

std::array<std::uint8_t, 44> createHeader(std::uint16_t numChannels, std::uint32_t sampleRate, std::uint16_t bitsPerSample, size_t samplesSize) {
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

    std::array<std::uint8_t, 44> bytes;

    auto writeHelper = [](const std::span<const std::uint8_t> &value, std::size_t offset, std::span<std::uint8_t> bytes) {
        if (offset + value.size() > bytes.size()) {
            throw std::out_of_range("writeHelper would write past end of buffer");
        }
        std::copy(value.begin(), value.end(), bytes.begin() + offset);
    };

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

double combine(const std::vector<Sound*> &sounds, size_t sampleIndex) {
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

        // TODO: this is wrong
        current *= amplitude->f(percent);

        combined += current;
    }

    combined /= sounds.size();

    return combined;
}

