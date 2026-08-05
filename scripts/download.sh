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
ln -s "$NUXSPACE" repos
mkdir -p arrangements/nuX
pushd arrangements/nuX
for thorn in "$NUXSPACE"/nuX_*; do
    ln -s "../../repos/$(basename "$thorn")" .
done
popd
