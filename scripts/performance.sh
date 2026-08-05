#!/bin/bash

set -euxo pipefail

export NUXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
cd "$WORKSPACE/Cactus"

export LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

SIMULATION_NAME=NuXPerformance
PARAMETER_FILE="$NUXSPACE/performance/nuX_M1_diffusion_cpu.par"

time ./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    create-run "$SIMULATION_NAME" --cores 2 --ppn-used 2 --num-threads 2 \
    --parfile "$PARAMETER_FILE"

PERFORMANCE_OUTPUT_DIR="$(./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    get-output-dir "$SIMULATION_NAME")"

if test -n "${GITHUB_ENV:-}"; then
    echo "PERFORMANCE_OUTPUT_DIR=$PERFORMANCE_OUTPUT_DIR" >>"$GITHUB_ENV"
fi

TIMER_REPORT="$(find "$PERFORMANCE_OUTPUT_DIR" -type f \
    -name 'TimerReport.000000.txt' -print -quit)"
ALL_TIMERS_TSV="$(find "$PERFORMANCE_OUTPUT_DIR" -type f \
    -name 'AllTimers.tsv' -print -quit)"
ALL_TIMERS_READABLE="$(find "$PERFORMANCE_OUTPUT_DIR" -type f \
    -name 'AllTimersReadable.txt' -print -quit)"

test -n "$TIMER_REPORT"
test -s "$TIMER_REPORT"
test -n "$ALL_TIMERS_TSV"
test -s "$ALL_TIMERS_TSV"
test -n "$ALL_TIMERS_READABLE"
test -s "$ALL_TIMERS_READABLE"

# These checks make this a stable instrumentation gate without using noisy
# wall-clock thresholds on a shared GitHub Actions runner.
grep -q '^Timer Report at iteration 16 ' "$TIMER_REPORT"
grep -q '^nuX_M1[[:space:]]*|' "$TIMER_REPORT"
grep -q '^ODESolvers[[:space:]]*|' "$TIMER_REPORT"
grep -q '^CarpetX[[:space:]]*|' "$TIMER_REPORT"
grep -q 'Total time for simulation' "$TIMER_REPORT"

TIMERS_TSV="$PERFORMANCE_OUTPUT_DIR/performance-timers.tsv"
SUMMARY_FILE="$PERFORMANCE_OUTPUT_DIR/performance-summary.md"

printf 'seconds\tthorn\troutine\n' >"$TIMERS_TSV"
awk -F '|' '
  $1 ~ /^(nuX_M1|ODESolvers|CarpetX)[[:space:]]*$/ {
    thorn = $1
    routine = $2
    seconds = $3
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", thorn)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", routine)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", seconds)
    if (seconds ~ /^[0-9]+([.][0-9]+)?$/)
      printf "%s\t%s\t%s\n", seconds, thorn, routine
  }
' "$TIMER_REPORT" | sort -t "$(printf '\t')" -k1,1nr >>"$TIMERS_TSV"

{
    printf '## nuX CPU performance smoke test\n\n'
    printf -- '- Case: M1 diffusion, two refinement levels, subcycling\n'
    printf -- '- Integrator: IMEX42L\n'
    printf -- '- Work: 16 iterations, one MPI rank, two CPU threads\n'
    printf -- '- Policy: timer integrity is gated; elapsed time is reported only\n\n'
    printf '```text\n'
    grep 'Total time for simulation' "$TIMER_REPORT" | tail -n 1
    printf '\nSlowest nuX/ODESolvers/CarpetX scheduled routines:\n'
    head -n 11 "$TIMERS_TSV"
    printf '```\n'
} >"$SUMMARY_FILE"

cat "$SUMMARY_FILE"
if test -n "${GITHUB_STEP_SUMMARY:-}"; then
    cat "$SUMMARY_FILE" >>"$GITHUB_STEP_SUMMARY"
fi
