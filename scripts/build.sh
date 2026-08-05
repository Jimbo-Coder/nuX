#!/bin/bash

set -euxo pipefail

export NUXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
cd "$WORKSPACE/Cactus"

ASTERX_SCRIPTS="$WORKSPACE/Cactus/arrangements/AsterX/scripts"
cp "$ASTERX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.cfg" \
    simfactory/mdb/optionlists/
cp "$ASTERX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.ini" \
    simfactory/mdb/machines/
cp "$ASTERX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.run" \
    simfactory/mdb/runscripts/
cp "$ASTERX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.sub" \
    simfactory/mdb/submitscripts/
cp "$ASTERX_SCRIPTS/defs.local.ini" simfactory/etc/
cp "$NUXSPACE/scripts/nux.th" .

sed -i \
    -e 's#/__w/AsterX/workspace#/__w/nuX/workspace#g' \
    -e 's#AsterX CI#nuX CI#g' \
    -e 's#github.com/jaykalinani/AsterX#github.com/jaykalinani/nuX#g' \
    "simfactory/mdb/machines/actions-$ACCELERATOR-$REAL_PRECISION.ini"
sed -i -e 's/thornlist       = asterx.th/thornlist       = nux.th/' \
    simfactory/etc/defs.local.ini

if command -v ccache >/dev/null 2>&1; then
    export CCACHE_DIR="${CCACHE_DIR:-$NUXSPACE/.ccache}"
    ccache --max-size=2G
    ccache --zero-stats || true
    sed -i -e 's/^CC = /CC = ccache /' -e 's/^CXX = /CXX = ccache /' \
        "simfactory/mdb/optionlists/actions-$ACCELERATOR-$REAL_PRECISION.cfg"
fi

time ./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    build -j "$(nproc)" sim 2>&1 | tee build.log

test -x exe/cactus_sim
command -v ccache >/dev/null 2>&1 && ccache --show-stats || true
