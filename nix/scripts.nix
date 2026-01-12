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
    
    echo "Running direct binary execution (no mpirun): $APP_DIR/bin/dc $@"
    exec $APP_DIR/bin/dc "$@"
  '';

  run-simgrid = pkgs.writeShellScriptBin "run-simgrid" ''
    BACKEND=$1
    SIZE=$2

    # Map backend to packages (simgrid ones)
    if [ "$BACKEND" == "cuda" ]; then
      APP_DIR=${packages.dc-simgrid-cuda}
      
      DRIVER_SANDBOX=$(mktemp -d)
      trap "rm -rf $DRIVER_SANDBOX" EXIT
      POSSIBLE_PATHS=("/usr/lib/x86_64-linux-gnu" "/usr/lib64" "/usr/lib/wsl/lib" "/usr/lib")
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
      if [ "$FOUND_DRIVER" -eq 1 ]; then export LD_LIBRARY_PATH="$DRIVER_SANDBOX:$LD_LIBRARY_PATH"; fi

    elif [ "$BACKEND" == "openmp" ]; then
      APP_DIR=${packages.dc-simgrid}
    else
      echo "Usage: run-simgrid [openmp|cuda] SIZE"
      exit 1
    fi

    absorption=2
    size_x=$(( SIZE - 2 * absorption - 8 ))
    size_y=$size_x
    size_z=$size_y
    
    # Defaults
    dx=1e-1; dy=1e-1; dz=1e-1; dt=1e-6; tmax=1e-4

    export OMP_NUM_THREADS=1

    ${pkgs.simgrid}/bin/smpirun -platform $APP_DIR/platform.xml \
      --cfg=smpi/display-timing:yes \
      --cfg=precision/timing:1e-9 \
      --cfg=tracing/precision:9 \
      --cfg=smpi/host-speed:auto \
      -trace --cfg=tracing/filename:dc.trace \
      -hostfile $APP_DIR/hostfile.txt \
      $APP_DIR/bin/dc \
      --size-x=$size_x --size-y=$size_y --size-z=$size_z \
      --absorption=$absorption \
      --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax \
      --output-file=./validation/predicted.dc

    ${pajeng}/bin/pj_dump -l 9 dc.trace | grep ^State > dc.csv
    ${rEnv}/bin/Rscript ./plot.R dc.csv
  '';

  run-dc-comparison = pkgs.writeShellScriptBin "run-dc-comparison" ''
    export PATH=${pkgs.cudatoolkit}/bin:$PATH
    size_x=52; size_y=52; size_z=52; absorption=2; dx=1e-1; dy=1e-1; dz=1e-1; dt=1e-6; tmax=1e-4
    
    # Use mpip variants for comparison (execution only)
    DC_OMP=${packages.dc-omp-mpip}
    DC_CUDA=${packages.dc-cuda-mpip}

    # Ensure MPIP is NOT active for ground truth if we reuse mpip binary, but cleaner to use plain if we had it.
    # But we reused mpip binary. Let's unset MPIP env var to be safe or set it to /dev/null if possible.
    # Or just ignore the output.
    export MPIP="-g" # Disable mpiP report generation

    echo "Running sequential version..."
    OMP_NUM_THREADS=1 ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none $DC_OMP/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/ground_truth.dc

    echo "Running OpenMP version..."
    ${pkgs.openmpi}/bin/mpirun -np 6 --bind-to none $DC_OMP/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/openmp_predicted.dc

    echo "Comparing..."
    mv ./validation/openmp_predicted.dc ./validation/predicted.dc
    Rscript ./validation/CompareResults.R 0

    # CUDA Setup
    DRIVER_SANDBOX=$(mktemp -d)
    trap "rm -rf $DRIVER_SANDBOX" EXIT
    POSSIBLE_PATHS=("/usr/lib/x86_64-linux-gnu" "/usr/lib64" "/usr/lib/wsl/lib" "/usr/lib")
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
    if [ "$FOUND_DRIVER" -eq 1 ]; then export LD_LIBRARY_PATH="$DRIVER_SANDBOX:$LD_LIBRARY_PATH"; fi

    echo "Running CUDA version..."
    ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none $DC_CUDA/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc

    echo "Comparing CUDA..."
    Rscript ./validation/CompareResults.R 1e-3
  '';

  tupi = pkgs.writeShellScriptBin "tupi" ''
    ${selectAppLogic}
    NP=$3
    MACHINEFILE=$4
    shift 4
    ARGS="$@"

    rm -f dc.mpiP* *.rst dc.trace dc.csv dc.output
    
    MPI_ARGS=""
    if [ "$PROFILE" == "mpip" ]; then 
       export MPIP="-k 2 -f ./dc.mpiP"
       MPI_ARGS="-x MPIP"
    fi

    echo "Running on Tupi: backend=$BACKEND profile=$PROFILE np=$NP machinefile=$MACHINEFILE"
    ${pkgs.openmpi}/bin/mpirun -np $NP \
      -machinefile $MACHINEFILE \
      --mca btl ^openib \
      --mca btl_tcp_if_include 192.168.0.30/24 \
      --bind-to none \
      $MPI_ARGS \
      $APP_DIR/bin/dc $ARGS | tee dc.output
    
    ${postProcessLogic}
  '';

  cei = pkgs.writeShellScriptBin "cei" ''
    ${selectAppLogic}
    NP=$3
    MACHINEFILE=$4
    shift 4
    ARGS="$@"

    rm -f dc.mpiP* *.rst dc.trace dc.csv dc.output
    
    MPI_ARGS=""
    if [ "$PROFILE" == "mpip" ]; then 
       export MPIP="-k 2 -f ./dc.mpiP"
       MPI_ARGS="-x MPIP"
    fi

    echo "Running on Cei: backend=$BACKEND profile=$PROFILE np=$NP machinefile=$MACHINEFILE"
    ${pkgs.openmpi}/bin/mpirun -np $NP \
      -machinefile $MACHINEFILE \
      --mca btl ^openib \
      --mca btl_tcp_if_include 192.168.0.30/24 \
      --bind-to none \
      $MPI_ARGS \
      $APP_DIR/bin/dc $ARGS | tee dc.output

    ${postProcessLogic}
  '';

  draco = pkgs.writeShellScriptBin "draco" ''
    ${selectAppLogic}
    NP=$3
    MACHINEFILE=$4
    shift 4
    ARGS="$@"

    rm -f dc.mpiP* *.rst dc.trace dc.csv dc.output
    
    MPI_ARGS=""
    if [ "$PROFILE" == "mpip" ]; then 
       export MPIP="-k 2 -f ./dc.mpiP"
       MPI_ARGS="-x MPIP"
    fi

    echo "Running on Draco: backend=$BACKEND profile=$PROFILE np=$NP machinefile=$MACHINEFILE"
    ${pkgs.openmpi}/bin/mpirun -np $NP \
      -machinefile $MACHINEFILE \
      --mca btl ^openib \
      --mca btl_tcp_if_include 192.168.0.30/24 \
      --bind-to none \
      $MPI_ARGS \
      $APP_DIR/bin/dc $ARGS | tee dc.output

    ${postProcessLogic}
  '';
}
