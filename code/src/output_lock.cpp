#include "precpack/output_lock.hpp"

#include <cerrno>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#endif

namespace precpack {

struct OutputLock::Impl {
    std::filesystem::path path;
#ifdef _WIN32
    HANDLE handle = INVALID_HANDLE_VALUE;
#else
    int descriptor = -1;
#endif
};

namespace {

[[nodiscard]] std::runtime_error lock_conflict(
    const std::filesystem::path& path) {
    return std::runtime_error(
        "another PrecPack process is using this output directory: " +
        path.parent_path().string());
}

}

OutputLock::OutputLock(const std::filesystem::path& output_directory)
    : impl_(std::make_unique<Impl>()) {
    std::filesystem::create_directories(output_directory);
    impl_->path = output_directory / ".precpack.lock";
#ifdef _WIN32
    impl_->handle = CreateFileW(
        impl_->path.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr,
        OPEN_ALWAYS, FILE_ATTRIBUTE_HIDDEN, nullptr);
    if (impl_->handle == INVALID_HANDLE_VALUE) {
        const DWORD code = GetLastError();
        if (code == ERROR_SHARING_VIOLATION || code == ERROR_LOCK_VIOLATION) {
            throw lock_conflict(impl_->path);
        }
        throw std::filesystem::filesystem_error(
            "cannot lock PrecPack output directory", impl_->path,
            std::error_code(static_cast<int>(code), std::system_category()));
    }
#else
    impl_->descriptor =
        open(impl_->path.c_str(), O_CREAT | O_RDWR | O_CLOEXEC, 0666);
    if (impl_->descriptor < 0) {
        throw std::filesystem::filesystem_error(
            "cannot open PrecPack output lock", impl_->path,
            std::error_code(errno, std::generic_category()));
    }
    if (flock(impl_->descriptor, LOCK_EX | LOCK_NB) != 0) {
        const int code = errno;
        close(impl_->descriptor);
        impl_->descriptor = -1;
        if (code == EWOULDBLOCK || code == EAGAIN) {
            throw lock_conflict(impl_->path);
        }
        throw std::filesystem::filesystem_error(
            "cannot lock PrecPack output directory", impl_->path,
            std::error_code(code, std::generic_category()));
    }
#endif
}

OutputLock::~OutputLock() {
#ifdef _WIN32
    if (impl_ && impl_->handle != INVALID_HANDLE_VALUE) {
        CloseHandle(impl_->handle);
    }
#else
    if (impl_ && impl_->descriptor >= 0) {
        static_cast<void>(flock(impl_->descriptor, LOCK_UN));
        static_cast<void>(close(impl_->descriptor));
    }
#endif
}

}
