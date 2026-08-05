#!/bin/bash

set -euxo pipefail

export NUXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$NUXSPACE/scripts/nux.th"

cd Cactus
mkdir -p arrangements/nuX
pushd arrangements/nuX
for thorn in "$NUXSPACE"/nuX_*; do
    test -f "$thorn/interface.ccl"
    ln -s "$thorn" .
    test -f "$(basename "$thorn")/interface.ccl"
done
popd
