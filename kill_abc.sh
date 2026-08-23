#!/bin/bash
# Kill everything abc_raster.R can spawn: the R/Rscript process itself, any
# future/furrr multisession workers (parallel:::.slaveRSOCK()), and any
# in-flight model/main simulation binaries launched via system2().
set -u

killed=0

kill_by_pattern() {
  local pattern="$1"
  local pids
  pids=$(pgrep -f "$pattern" || true)
  if [ -n "$pids" ]; then
    echo "Killing ($pattern): $pids"
    kill -9 $pids
    killed=1
  fi
}

kill_by_pattern "R/abc_raster.R"
kill_by_pattern "parallel:::.slaveRSOCK"
kill_by_pattern "model/main"
kill_by_pattern "\./main"

if [ "$killed" -eq 0 ]; then
  echo "Nothing matching abc_raster.R, future workers, or main found."
else
  echo "Done."
fi
