{ pkgs, packages, pajeng, rEnv, mpiP, akypuera }:
let 
  # Helper for backend selection and CUDA env setup
  selectAppLogic = ''
    BACKEND=$1
    PROFILE=$2
    
    # Defaults
    if [ -z "$BACKEND" ]; then BACKEND="openmp"; fi
    if [ -z "$PROFILE" ]; then PROFILE="mpip"; fi

    APP_DIR=""
    if [ "$BACKEND" == "openmp" ]; then
      if [ "$PROFILE" == "mpip" ]; then APP_DIR=${packages.dc-omp-mpip}; fi
      if [ "$PROFILE" == "akypuera" ]; then APP_DIR=${packages.dc-omp-aky}; fi
    elif [ "$BACKEND" == "cuda" ]; then
       # CUDA Driver Sandbox
       DRIVER_SANDBOX=$(mktemp -d)
       trap "rm -rf $DRIVER_SANDBOX" EXIT
       POSSIBLE_PATHS=(
         "/usr/lib/x86_64-linux-gnu"
         "/usr/lib64"
         "/usr/lib/wsl/lib"
         "/usr/lib"
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

      if [ "$PROFILE" == "mpip" ]; then APP_DIR=${packages.dc-cuda-mpip}; fi
      if [ "$PROFILE" == "akypuera" ]; then APP_DIR=${packages.dc-cuda-aky}; fi
    fi
    
    if [ -z "$APP_DIR" ]; then
      echo "Error: Invalid backend ($BACKEND) or profile ($PROFILE)"
      echo "Supported backends: openmp, cuda"
      echo "Supported profiles: mpip, akypuera"
      exit 1
    fi
  '';

  # Helper for post-processing logic
  postProcessLogic = ''
    echo "Starting post-processing for profile: $PROFILE"
    if [ "$PROFILE" == "mpip" ]; then
       MPIP_FILE=$(ls *.mpiP* 2>/dev/null | head -n 1)
       if [ -n "$MPIP_FILE" ]; then
         echo "Processing mpiP file: $MPIP_FILE"
         ${rEnv}/bin/Rscript ./plot.R $MPIP_FILE
       elif [ -f "dc.output" ]; then
         echo "No .mpiP file found. Attempting to parse dc.output..."
         ${rEnv}/bin/Rscript ./plot.R dc.output
       else
         echo "Warning: No mpiP output found (checked *.mpiP* and dc.output)."
         ls -la
       fi
    elif [ "$PROFILE" == "akypuera" ]; then
       if ls rastro-*.rst 1> /dev/null 2>&1; then
         echo "Converting Akypuera traces..."
         ${akypuera}/bin/aky_converter > dc.trace
         ${pajeng}/bin/pj_dump -l 9 dc.trace | grep ^State > dc.csv
         ${rEnv}/bin/Rscript ./plot.R dc.csv
       else
         echo "Warning: No Akypuera traces found (rastro-*.rst)."
       fi
    fi
  '';

in
{
  run-dc = pkgs.writeShellScriptBin "run-dc" ''
    ${selectAppLogic}
    shift 2
    ARGS="$@"
    
    # Cleanup
    rm -f dc.mpiP* *.rst dc.trace dc.csv dc.output

    # Environment
    MPI_ARGS=""
    if [ "$PROFILE" == "mpip" ]; then
       export MPIP="-k 2 -f ./dc.mpiP"
       MPI_ARGS="-x MPIP"
    fi

    echo "Running $APP_DIR/bin/dc with args: $ARGS"
    ${pkgs.openmpi}/bin/mpirun -np 6 --bind-to none $MPI_ARGS $APP_DIR/bin/dc $ARGS | tee dc.output

    ${postProcessLogic}
  '';

  dc = pkgs.writeShellScriptBin "dc" ''
    ${selectAppLogic}
    shift 2
    
    if [ "$PROFILE" == "mpip" ]; then
       export MPIP="-k 2 -f $(pwd)/dc.mpiP"
    fi

    echo "Running direct binary execution (no mpirun): $APP_DIR/bin/dc $@"
    exec $APP_DIR/bin/dc "$@"
  '';

  run-simgrid = pkgs.writeShellScriptBin "run-simgrid" ''
    BACKEND=$1
    shift 1
    ARGS="$@"

    APP_DIR=""
    if [ "$BACKEND" == "openmp" ]; then
      APP_DIR=${packages.dc-simgrid}
    elif [ "$BACKEND" == "cuda" ]; then
      APP_DIR=${packages.dc-simgrid-cuda}
      
      DRIVER_SANDBOX=$(mktemp -d)
      trap "rm -rf $DRIVER_SANDBOX" EXIT
      POSSIBLE_PATHS=(
        "/usr/lib/x86_64-linux-gnu"
        "/usr/lib64"
        "/usr/lib/wsl/lib"
        "/usr/lib"
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
        echo "WARNING: Could not find host libcuda.so.1. Simulation might crash."
      fi
    else
      echo "Usage: run-simgrid [openmp|cuda] [args...]"
      exit 1
    fi

    rm -f dc.trace dc.csv smpi.html smpi.png

    export OMP_NUM_THREADS=1

    echo "Running SimGrid ($BACKEND) with args: $ARGS"
    ${pkgs.simgrid}/bin/smpirun -platform $APP_DIR/platform.xml \
      --cfg=smpi/display-timing:yes \
      --cfg=precision/timing:1e-9 \
      --cfg=tracing/precision:9 \
      --cfg=smpi/host-speed:auto \
      -trace --cfg=tracing/filename:dc.trace \
      -hostfile $APP_DIR/hostfile.txt \
      $APP_DIR/bin/dc $ARGS

    if [ -f "dc.trace" ]; then
        ${pajeng}/bin/pj_dump -l 9 dc.trace | grep ^State > dc.csv
        ${rEnv}/bin/Rscript ./plot.R dc.csv
    else
        echo "No trace generated."
    fi
  '';

  run-dc-comparison = pkgs.writeShellScriptBin "run-dc-comparison" ''
