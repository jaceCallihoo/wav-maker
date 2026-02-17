#include <arpa/inet.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <array>

// struct {
//     
// } wavData;

// void writeData(wavData input) std::vector<std::uint8_t> {
// 
// }


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


class Wav {
    // ChunkID          constant
    // ChunkSize        derived
    // Format           constant
    // Subchunk1ID      constant
    // Subchunk1Size    constant maybe?
    // AudioFormat      constant (1 for PCM)
    // NumChannels      variable (1 = mono, 2 = stereo)
    // SampleRate       variable (8000, 44100, etc)
    // ByteRate         derived
    // BlockAlign       derived
    // BitsPerSample    variable (8 = 8 bits, 16 = 16 bits)
    // ExtraParmasSize  other (if PCM, then doesn't exist)
    // ExtraParams      other (space for extra parameters)
    // Subchunk2ID      constant ("data")
    // Subchunk2Size    derived (NumSamples * NumChannels * BitsPerSample/8)
    // Data             variable

    static constexpr std::array<std::uint8_t, 4> chunkID = {'R', 'I', 'F', 'F'};
    // chunkSize
    static constexpr std::array<std::uint8_t, 4> format = {'W', 'A', 'V', 'E'};
    static constexpr std::array<std::uint8_t, 4> subchunk1ID = {'f', 'm', 't', ' '};
    static constexpr std::array<std::uint8_t, 4> subchunk1Size = toLittleEndian<std::uint32_t>(16);
    static constexpr std::array<std::uint8_t, 4> audioFormat = toLittleEndian<std::uint32_t>(1);
    int numChannels;
    int sampleRate;
    // byteRate
    // blockAlign
    int bitsPerSample;
    // extraParamsSize
    // extraParams
    static constexpr std::array<std::uint8_t, 4> subchunk2ID = {'d', 'a', 't', 'a'};
    // subchunk2Size
    std::vector<std::uint8_t> data;

    int deriveChunkSize() {
        return 36 + deriveSubChunk2Size();
    }

    int deriveByteRate() {
        return sampleRate * numChannels * bitsPerSample / 8;
    }

    int deriveBlockAlign() {
        return numChannels * bitsPerSample / 8;
    }

    int deriveSubChunk2Size() {
        return data.size();
    }
public:
    std::vector<std::uint8_t> toBytes() {
        return std::vector<std::uint8_t>();
    }
};

int main() {
    std::cout << "hello world" << std::endl;

    //////////////////////////////////////////////////////////////////////
    
    std::vector<std::uint8_t> buffer;
    buffer.insert(buffer.end(), {'R','I','F','F'});

    std::ofstream file("output.txt", std::ios::binary);
    if (!file) {
        std::cerr << "failed to open file\n";
        return 1;
    }

    file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    // optional, this is called when the file leaves scope
    file.close();

    return 0;
}
