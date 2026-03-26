{
  pkgs,
  packages,
  pajeng,
  rEnv,
  akypuera,
}: let
  # Helper for backend selection and CUDA env setup
  selectAppLogic = ''
    BACKEND=$1

    # Default
    if [ -z "$BACKEND" ]; then BACKEND="openmp"; fi

    APP_DIR=""
    if [ "$BACKEND" == "openmp" ]; then
      APP_DIR=${packages.dc-omp-aky}
    elif [ "$BACKEND" == "cuda" ]; then
       # CUDA Driver Sandbox
       DRIVER_SANDBOX=$(mktemp -d)
       trap "rm -rf $DRIVER_SANDBOX" EXIT
       POSSIBLE_PATHS=(
         "/usr/lib/x86_64-linux-gnu"
         "/usr/lib64"
         "/usr/lib/wsl/lib"
         "/usr/lib"
         "/run/opengl-driver/lib"
       )
       FOUND_DRIVER=0
       for libdir in "''${POSSIBLE_PATHS[@]}"; do
         if [ -e "$libdir/libcuda.so.1" ]; then
           echo "Found host CUDA driver in: $libdir"
           ln -sf "$libdir/libcuda.so.1" "$DRIVER_SANDBOX/libcuda.so.1"
           ln -sf "$libdir/libcuda.so.1" "$DRIVER_SANDBOX/libcuda.so"
           if [ -e "$libdir/libnvidia-ptxjitcompiler.so.1" ]; then
               ln -sf "$libdir/libnvidia-ptxjitcompiler.so.1" "$DRIVER_SANDBOX/libnvidia-ptxjitcompiler.so.1"
           fi
           FOUND_DRIVER=1
           break
         fi
       done
       if [ "$FOUND_DRIVER" -eq 1 ]; then
         export LD_LIBRARY_PATH="$DRIVER_SANDBOX:$LD_LIBRARY_PATH"
       else
         echo "WARNING: Could not find host libcuda.so.1."
       fi

       APP_DIR=${packages.dc-cuda-aky}
    fi

    if [ -z "$APP_DIR" ]; then
      echo "Error: Invalid backend ($BACKEND)"
      echo "Supported backends: openmp, cuda"
      exit 1
    fi
  '';

  # Helper for post-processing logic
  postProcessLogic = ''
    echo "Starting post-processing..."
    if ls rastro-*.rst 1> /dev/null 2>&1; then
      echo "Converting Akypuera traces..."
      ${akypuera}/bin/aky_converter *.rst > dc.trace
      ${pajeng}/bin/pj_dump -z -l 9 dc.trace | grep ^State > dc.csv
      ${rEnv}/bin/Rscript ./plot.R dc.csv
    else
      echo "Warning: No Akypuera traces found (rastro-*.rst)."
    fi
  '';

  runSimgridPlatformCuda = pkgs.writeShellScriptBin "run-simgrid-platform-cuda" ''
    export PATH=${pkgs.lib.makeBinPath [
      pkgs.coreutils
      pkgs.gnugrep
      pkgs.simgrid
      packages.dc-simgrid-platform
      pajeng
      rEnv
    ]}:$PATH

    if [ "$#" -lt 6 ]; then
      echo "Usage: run-simgrid-platform-cuda NUM_HOSTS NET_BW NET_LAT GPU_BW GPU_LAT GPU_POWER [PROGRAM_ARGS...]"
      exit 1
    fi

    NUM_HOSTS=$1
    NET_BW=$2
    NET_LAT=$3
    shift 3
    ARGS="$@"

    export PLATFORM_NUM_HOSTS=$NUM_HOSTS
    export PLATFORM_NET_BW=$NET_BW
    export PLATFORM_NET_LAT=$NET_LAT
    export PLATFORM_HOSTFILE=simgrid-config/hostfile.txt

    # CUDA Setup
    DRIVER_SANDBOX=$(mktemp -d)
    trap "rm -rf $DRIVER_SANDBOX" EXIT
    POSSIBLE_PATHS=(
      "/usr/lib/x86_64-linux-gnu"
      "/usr/lib64"
      "/usr/lib/wsl/lib"
      "/usr/lib"
      "/run/opengl-driver/lib"
    )
    FOUND_DRIVER=0
    for libdir in "''${POSSIBLE_PATHS[@]}"; do
      if [ -e "$libdir/libcuda.so.1" ]; then
        ln -sf "$libdir/libcuda.so.1" "$DRIVER_SANDBOX/libcuda.so.1"
        ln -sf "$libdir/libcuda.so.1" "$DRIVER_SANDBOX/libcuda.so"
        if [ -e "$libdir/libnvidia-ptxjitcompiler.so.1" ]; then
            ln -sf "$libdir/libnvidia-ptxjitcompiler.so.1" "$DRIVER_SANDBOX/libnvidia-ptxjitcompiler.so.1"
        fi
        FOUND_DRIVER=1
        break
      fi
    done
    if [ "$FOUND_DRIVER" -eq 1 ]; then
      export LD_LIBRARY_PATH="$DRIVER_SANDBOX:$LD_LIBRARY_PATH"
    fi

    # Generate Hostfile (using pre-built tool from PATH)
    generate_artifacts > /dev/null

    # Run Simulation
    export OMP_NUM_THREADS=1
    rm -f dc.trace dc.csv smpi.html smpi.png

    # Use smpirun from PATH (pkgs.simgrid)
    # Use libplatform.so from packages.dc-simgrid-platform (via PATH/../lib implicit or explicit path)

    PLATFORM_LIB=${packages.dc-simgrid-platform}/lib/libplatform.so
    DC_BIN=${packages.dc-simgrid-cuda}/bin/dc

    ${pkgs.simgrid}/bin/smpirun \
      -platform $PLATFORM_LIB \
      -hostfile simgrid-config/hostfile.txt \
      -np $NUM_HOSTS \
      --cfg=smpi/display-timing:yes \
      --cfg=precision/timing:1e-9 \
      --cfg=tracing/precision:9 \
      --cfg=smpi/shared-malloc:global \
      --cfg=smpi/host-speed:"1f" \
      --cfg=smpi/display-allocs:yes \
      --cfg=smpi/auto-shared-malloc-thresh:1073741824 \
      -trace --cfg=tracing/filename:dc.trace \
      $DC_BIN $ARGS 2>&1 | tee sim.log

    # Generate Stats
    if [ -f "dc.trace" ]; then
        pj_dump -l 9 dc.trace | grep ^State > dc.csv
        ${rEnv}/bin/Rscript ./get_metrics.R
    fi
  '';
  runSimGridExperiments = ''
    OUTPUT_CSV="simulation_results.csv"
    rm -rf "$OUTPUT_CSV" "sim_traces" "plots"
    echo "run,problem_size,mpi_time,computation_time,total_time" > $OUTPUT_CSV

    ANALYSIS_DIR=$1
    PROBLEM_SIZES=$2

    IFS=',' read -r -a PROBLEM_SIZES <<< "$PROBLEM_SIZES"

    echo "Running experiments (sizes $PROBLEM_SIZES)..."

    for size in "''${PROBLEM_SIZES[@]}"; do
      echo "  Problem Size: $size"
      absorption=2
      cuda_size=$(( size - 8 - 2 * absorption ))

      output=$(${runSimgridPlatformCuda}/bin/run-simgrid-platform-cuda $NUM_HOSTS $NET_BW $NET_LAT  \
               --size-x=$cuda_size --size-y=$cuda_size --size-z=$cuda_size --absorption=$absorption --dx=1e-1 --dy=1e-1 --dz=1e-1 --dt=1e-6 --time-max=1e-3 --output-file=./validation/predicted.dc)

      total_time=$(echo "$output" | grep "Total time:" | awk '{print $4}' | tr -d '"')
      mpi_time=$(echo "$output" | grep "MPI time:" | awk '{print $4}' | tr -d '"')
      comp_time=$(echo "$output" | grep "Computation time:" | awk '{print $4}' | tr -d '"')

      if [ -n "$total_time" ]; then
        echo "$i,$size,$mpi_time,$comp_time,$total_time" >> $OUTPUT_CSV
        mkdir -p "sim_traces/$size"
        cp dc.csv "sim_traces/$size/dc.csv"
      else
        echo "    Failed to parse metrics for run size $size:"
        echo "    $output"
      fi
    done

    echo "Comparison vs Real Data:"
    if [ -d "$ANALYSIS_DIR" ]; then
      ${rEnv}/bin/Rscript ./validation/compare_sim_real.R $ANALYSIS_DIR sim_traces
    else
      echo "Warning: Analysis directory not found: $ANALYSIS_DIR"
      echo "Skipping comparison with real data."
    fi
  '';

  runPotiExperiments = pkgs.writeShellScriptBin "run-poti-experiments" ''
    export PATH=${pkgs.lib.makeBinPath [
      pkgs.coreutils
      pkgs.gawk
      pkgs.gnugrep
      rEnv
      pkgs.pandoc
      akypuera
      pajeng
    ]}:$PATH

    NUM_HOSTS=5
    NET_BW="937Mbps"
    NET_LAT="22.7us"

    ${runSimGridExperiments}
  '';
in {
  run-dc = pkgs.writeShellScriptBin "run-dc" ''
    ${selectAppLogic}
    np=$2
    shift 2
    ARGS="$@"

    # Cleanup
    rm -f *.rst dc.trace dc.csv dc.output

    echo "Running $APP_DIR/bin/dc with args: $ARGS"
    ${pkgs.openmpi}/bin/mpirun -np $np --bind-to none $APP_DIR/bin/dc $ARGS | tee dc.output

    ${postProcessLogic}
  '';

  dc = pkgs.writeShellScriptBin "dc" ''
    ${selectAppLogic}
    shift 1

    echo "Running direct binary execution (no mpirun): $APP_DIR/bin/dc $@"
    exec $APP_DIR/bin/dc "$@"
  '';

  run-simgrid-platform-cuda = runSimgridPlatformCuda;
  run-poti-experiments = runPotiExperiments;
}
