#!/bin/bash
# Temporary benchmark script — delete after use

CP=0.8   # commuting proportion — change this to rerun with a different value

SCALES=(Raster220 Raster220_10k Raster220_20k Raster220_100k Raster660 Raster660_10k Raster660_20k Raster660_100k)
FITTED=(Raster220 Raster220     Raster220      Raster220      Raster660 Raster660     Raster660     Raster660)
T2=1
REPS=2
RESULTS="../output/benchmark_times.csv"

# ── One-time migration: rename old-style files (no cp in name) ────────────────
for f in ../output/benchmark_*.csv; do
    base=$(basename "$f")
    if [[ "$base" != benchmark_cp* && "$base" != benchmark_times.csv ]]; then
        mv "$f" "../output/benchmark_cp0.5_${base#benchmark_}"
        echo "Renamed $base -> benchmark_cp0.5_${base#benchmark_}"
    fi
done
for f in ../output/benchmark_*.netfil; do
    base=$(basename "$f")
    if [[ "$base" != benchmark_cp* ]]; then
        mv "$f" "../output/benchmark_cp0.5_${base#benchmark_}"
    fi
done
# Update times CSV to add cp column if still in old format
if [ -f "$RESULTS" ] && [ "$(head -1 "$RESULTS")" = "scale,rep,seconds" ]; then
    tmpfile=$(mktemp)
    echo "cp,scale,rep,seconds" > "$tmpfile"
    tail -n +2 "$RESULTS" | sed 's/^/0.5,/' >> "$tmpfile"
    mv "$tmpfile" "$RESULTS"
    echo "Migrated benchmark_times.csv to include cp column"
fi

# ── Patch params.h and recompile ─────────────────────────────────────────────
sed -i '' -E "s/(constexpr double COMMUTING_PROP[[:space:]]*=[[:space:]]*)[0-9.]+/\1${CP}/" params.h
echo "Compiling with COMMUTING_PROP = ${CP}"
g++ -std=c++20 -O2 -o main \
    agent.cpp dynamic_network.cpp initialise_pop.cpp main.cpp mda.cpp \
    mobility.cpp rand_func.cpp rng.cpp sim.cpp statistics.cpp write_netfil_log.cpp
echo "Done compiling"

mkdir -p ../output
if [ ! -f "$RESULTS" ]; then
    echo "cp,scale,rep,seconds" > "$RESULTS"
fi

for i in "${!SCALES[@]}"; do
    scale="${SCALES[$i]}"
    fitted="${FITTED[$i]}"

    for rep in $(seq 1 $REPS); do
        echo "=== cp=${CP} $scale rep $rep ==="

        ./clean_inputs.sh
        cp ../data/Scales/$scale/* ../data/
        cp ../data/Fitted/$fitted/Theta_$T2/tran_params.csv ../data/
        echo "$scale" > ../data/current_scale.txt

        start=$(python3 -c "import time; print(time.time())")
        ./main "benchmark_cp${CP}_${scale}_rep${rep}.csv"
        end=$(python3 -c "import time; print(time.time())")
        elapsed=$(python3 -c "print(round($end - $start, 2))")

        echo "${CP},$scale,$rep,$elapsed" >> "$RESULTS"
        echo "  -> ${elapsed}s"
    done
done

echo "Done. Results in $RESULTS"
echo "Benchmark outputs to clean up: ../output/benchmark_cp*.csv ../output/benchmark_cp*.netfil"
