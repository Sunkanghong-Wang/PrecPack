#pragma once

#include <filesystem>
#include <memory>

namespace precpack {

class OutputLock {
public:
    explicit OutputLock(const std::filesystem::path& output_directory);
    ~OutputLock();

    OutputLock(const OutputLock&) = delete;
    OutputLock& operator=(const OutputLock&) = delete;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}
