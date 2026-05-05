REM Get vcpkg
git clone --branch 2026.02.27 git@github.com:microsoft/vcpkg.git
set VCPKG_ROOT=%CD%\vcpkg
bootstrap-vcpkg.bat -S . -B build -I install -Y