#pragma once

#include <cstdlib>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace precpack::internal {

[[nodiscard]] inline std::optional<std::string> environment_value(
    const char* name) {
#if defined(_WIN32)
    char* raw_value = nullptr;
    std::size_t length = 0U;
    const errno_t error = _dupenv_s(&raw_value, &length, name);
    std::unique_ptr<char, decltype(&std::free)> value(
        raw_value, &std::free);
    if (error != 0) {
        throw std::runtime_error(
            "cannot read environment variable: " + std::string(name));
    }
    if (!value || length <= 1U) {
        return std::nullopt;
    }
    return std::string(value.get());
#else
    const char* value = std::getenv(name);
    if (value == nullptr || *value == '\0') {
        return std::nullopt;
    }
    return std::string(value);
#endif
}

}
