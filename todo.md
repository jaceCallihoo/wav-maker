

setup testing framework
- gtest
setup linting and formatting
setup a package manager?
add utilities
- static analysis tool
    - cpp check
    - clang static analyzer
- code coverage tools
- docker
- dev container?

make sure I'm using c++ 20
choose a compiler
- what are the tradeoffs with the options?

plan on how to do depenedncy management
- static vs dynamic linking


Need to create a data structure which represents the wav file format that can have it's fields written to arbitrarily
Options for the data structure
### just a vec\<uint8_t\>
would be of size 44 + sizeof(data size)

### a class
has a to bytes function which can be written to file
can apply bytes in order, but maybe this is bad and the offset should always be targeted when writing
I like the class solution better
- because it's prabably annoying to manage the buffer size
- should maybe try the other solution later



###


define a bunch of different types of waves
- sin wave
- square wave

find out how to combine different types of waves

create a json parser
```json
{
    "durationMs": 5000,     // optional, can be derived from the waves delays and duration, can overwrite the duration of the other fields maybe?
    "sampleRate": 44100,    // required (all waves should be the same for now) also maybe should have a default
    "bitsPerSample": 16,    // maybe this should have a default later
    "channels": 1,          // default 1
    "waves": [
        {
            "type": "square",       // required
            "frequency": 100,       // required
            "amplitude": 16000      // required
            "delayMs": 0,           // default: 0
            "durationMs": 3000      // required
        },
        {
            "type": "sin",
            "frequency": 100,
            "amplitude": 16000
            "delayMs": 0,
            "durationMs": 3000
        },
    ]
}
```


create a wav file visualizer



    // uint8_t streakValue1 = buffer[0];
    // uint8_t streakValue2 = buffer[1];
    // uint8_t streakCount = 1;
    // for (int i = 0; i < buffer.size() + 2; i += 2) {
    //     if (buffer[i] != streakValue1 || buffer[i+1] != streakValue2) {
    //         std::cout << std::format("{}: {:x} {:x}\n", streakCount, int(streakValue1), int(streakValue2));
    //         streakValue1 = buffer[i];
    //         streakValue2 = buffer[i+1];
    //         streakCount = 1;
    //     } else {
    //         streakCount++;
    //     }
    // }

