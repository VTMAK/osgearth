/* osgEarth
* Copyright 2025 Pelican Mapping
* MIT License
*/

#define CATCH_CONFIG_RUNNER
#include <osgEarth/catch.hpp>
#include <osgEarth/FileUtils>
#include <osgDB/ConvertUTF>
#include <string>
#include <vector>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#elif defined(__linux__)
#  include <unistd.h>
#elif defined(__APPLE__)
#  include <mach-o/dyld.h>
#endif

namespace osgEarth { namespace Tests
{
    std::string executablePath;
} }

int main(int argc, char* argv[])
{
#ifdef _WIN32
    std::vector<wchar_t> buffer(32768u, L'\0');
    const DWORD length = GetModuleFileNameW(nullptr, buffer.data(),
        static_cast<DWORD>(buffer.size()));
    if (length > 0u && length < buffer.size())
        osgEarth::Tests::executablePath = osgDB::convertUTF16toUTF8(
            buffer.data(), length);
#elif defined(__linux__)
    std::vector<char> buffer(4096u, '\0');
    const ssize_t length = ::readlink("/proc/self/exe", buffer.data(), buffer.size() - 1u);
    if (length > 0)
        osgEarth::Tests::executablePath.assign(
            buffer.data(), static_cast<std::size_t>(length));
#elif defined(__APPLE__)
    uint32_t size = 0u;
    if (_NSGetExecutablePath(nullptr, &size) == -1)
    {
        std::vector<char> buffer(size, '\0');
        if (_NSGetExecutablePath(buffer.data(), &size) == 0)
            osgEarth::Tests::executablePath = buffer.data();
    }
#endif

    if (osgEarth::Tests::executablePath.empty())
        osgEarth::Tests::executablePath = osgEarth::Util::getAbsolutePath(argv[0]);

    return Catch::Session().run(argc, argv);
}
