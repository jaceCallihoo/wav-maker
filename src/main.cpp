#include <arpa/inet.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <array>
#include <span>

// TODO: shoulde the function names be camel case? Check google formatting standards
// TODO: check the class field names as well

template <typename T>
constexpr std::array<uint8_t, sizeof(T)> toLittleEndian(T value) {
    std::array<uint8_t, sizeof(T)> out{};
    for (size_t i = 0; i < sizeof(T); ++i) {
        out[i] = static_cast<char>((value >> (8 * i)) & 0xFF);
    }
    return out;
}

void writeHelper(std::span<const std::uint8_t> value, std::size_t offset, std::vector<std::uint8_t>& bytes) {
    if (offset + value.size() > bytes.size()) {
        throw std::out_of_range("writeHelper would write past end of buffer");
    }

    std::copy(value.begin(), value.end(), bytes.begin() + offset);
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

// sampleRate: how many samples of the wave are there in a seccond
// bitsPerSample: how many bits are used to record a sample
//


template <typename T>
std::vector<std::uint8_t> allToLittleEndian(std::span<T> values) {
    std::vector<std::uint8_t> bytes(values.size() * sizeof(T));

    for (std::size_t i = 0; i < values.size(); i++) {
        size_t offset = i * sizeof(T);
        writeHelper(values[i], offset, bytes);
    }

    return bytes;
}

// maybe this should return a std::int16_t instead of a std::uint8_t, it should be converted to little endian later
std::vector<std::uint8_t> squareWave(int durationMs, int sampleRate, int bitsPerSample, int frequency, std::int16_t high, std::int16_t low) {
    // assert frequency > 0
    // assert cycleDuration > 0
    // assert duration >= 0
    // assert sampleRate >= 0
    // assert bitsPerSample >= 8 && bitsPerSample % 8 == 0
    // assert high <= 2^(bitsPerSample - 1) - 1 && low >= -2^(bitsPerSample - 1)
    // assert bistPerSample / 8 == sizeof(high) && bitsPerSample / 8 == sizeof(low)

    const int bytesPerSample = bitsPerSample / 8;
    const int numSamples = (durationMs * sampleRate) / 1000;
    const int cycleDuration = sampleRate / frequency;
    int cycleCount = 0;
    bool isHigh = true;

    std::vector<std::uint8_t> sound(numSamples * bytesPerSample);

    for (std::size_t i = 0; i < numSamples; i++) {
        if (!(cycleCount < cycleDuration)) {
            isHigh = !isHigh;
            cycleCount = 0;
        }
        cycleCount++;

        const std::size_t offset = i * bytesPerSample;
        writeHelper(toLittleEndian(isHigh ? high : low), offset, sound);
    }

    return sound;
}

// need to take in the sample rate, or do I? yes becaues theyre signed i think 
std::vector<std::uint8_t> addWaves(std::vector<std::uint8_t> a, std::vector<std::uint8_t> b) {
    // assert a.size() == b.size()


}

std::vector<std::uint8_t> createSound(int numChannels, int sampleRate, int bitsPerSample) {
    const int durationSeconds = 1;
    const std::int16_t up = 16383;
    const std::int16_t down = -16383;

    std::vector<std::uint8_t> sound(durationSeconds * sampleRate * bitsPerSample / 8);

    bool isUp = true;
    // bool mod = true;
    for (int i = 0; i < durationSeconds * sampleRate; i++) {
        if (i % 100 == 0) {
            isUp = !isUp;
        }

        std::int16_t sample;
        if (isUp) {
            sample = up;
        } else {
            sample = down;
        }

        /*
        if (i % 300 == 0) {
            mod = !mod;
        }

        if (mod) {
            sample = sample / 2;
        }
        */

        auto serializedSample = toLittleEndian(sample);
        const size_t offset = i * bitsPerSample / 8;
        writeHelper(serializedSample, offset, sound);
    }

    return sound;
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
    std::vector<std::uint8_t> data = squareWave(2000, sampleRate, bitsPerSample, 440, 16383, -16383);

    Wav wave(numChannels, sampleRate, bitsPerSample, data);


    std::vector<std::uint8_t> buffer = wave.toBytes();


    //////////////////////////////////////////////////////////////////////
    // write to file

    // std::vector<std::uint8_t> buffer;
    // buffer.insert(buffer.end(), {'R','I','F','F'});

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
