#include <cstdint>
#include <vector>
#include <fstream>

//////////////// ///////////// ///////////// /////////////

#include <span>
#include <cstdint>

template<typename Sink>
concept ByteSink = requires(Sink& s, std::span<const std::uint8_t> data) {
    { s.write(data) } -> std::same_as<void>;
};

template<ByteSink Sink>
void writeData(Sink& sink)
{
    std::array<std::uint8_t, 2> a{0x01, 0x02};
    std::array<std::uint8_t, 4> b{0x03, 0x04, 0x05, 0x06};

    sink.write(a);
    sink.write(b);
}

#include <ostream>

struct OstreamSink {
    std::ostream& out;

    void write(std::span<const std::uint8_t> data)
    {
        out.write(
            reinterpret_cast<const char*>(data.data()),
            data.size()
        );
    }
};

//////////////////////////////////////////////////////////////////////////////

struct IBuffer {
    virtual void write(uint8_t byte) = 0;
    virtual void write(const uint8_t* data, size_t size) {
        for (size_t i = 0; i < size; ++i)
            write(data[i]);
    }
    virtual ~IBuffer() = default;
};

struct FileBuffer : IBuffer {
    std::ofstream out;
    FileBuffer(const std::string& filename)
        : out(filename, std::ios::binary) {}
    
    void write(uint8_t byte) override {
        out.put(static_cast<char>(byte));
    }
};

struct MemoryBuffer : IBuffer {
    std::vector<uint8_t> data;
    
    void write(uint8_t byte) override {
        data.push_back(byte);
    }
};
