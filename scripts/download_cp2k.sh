# Bash script for downloading CP2K

set -euo pipefail

ARCH="$1"
VERSION="$2"

# After the 9.1 release CP2K switched to a `<year>.<suffix>` version scheme (e.g. `2021.1`)
if [[ $VERSION =~ [0-9][0-9][0-9][0-9]\.[0-9]+ ]]; then
    PLAT="Linux-gnu"
    VERSION_LONG=v"$VERSION"
else
    PLAT="Linux"
    VERSION_LONG=v"$VERSION".0
fi

echo ::group::"Get CP2K $VERSION data"
git clone https://github.com/cp2k/cp2k -q
cd cp2k
git -c advice.detachedHead=false checkout tags/$VERSION_LONG
echo -e ::endgroup::\n

echo ::group::"Installing CP2K $ARCH $VERSION binaries"
curl -Lsf https://github.com/cp2k/cp2k/releases/download/$VERSION_LONG/cp2k-$VERSION-$PLAT-$ARCH.ssmp -o /usr/local/bin/cp2k.ssmp
chmod u+rx /usr/local/bin/cp2k.ssmp
echo -e ::endgroup::\n
