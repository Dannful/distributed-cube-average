{
  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.11";
    nixpkgs-cuda.url = "github:nixos/nixpkgs/nixos-24.11";
  };
  outputs = {
    self,
    nixpkgs,
    nixpkgs-cuda,
    utils,
  }:
    utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {
        inherit system;
        config = {
          allowUnfree = true;
        };
      };
      pkgs-cuda = import nixpkgs-cuda {
        inherit system;
        config = {
          allowUnfree = true;
        };
      };

      akypuera = import ./akypuera.nix {inherit pkgs;};

      pajeng = import ./pajeng.nix {inherit pkgs;};

      dc = pkgs.stdenv.mkDerivation {
        pname = "dc";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = with pkgs; [gnumake openmpi akypuera];
        buildPhase = "make all";
        installPhase = ''
          mkdir -p $out/bin
          cp bin/dc $out/bin
        '';
      };
      dc-simgrid = pkgs.stdenv.mkDerivation {
        pname = "dc";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = with pkgs; [gnumake openmpi simgrid];
        buildPhase = "make all BACKEND=simgrid";
        installPhase = ''
          mkdir -p $out/bin
          cp bin/dc $out/bin
          cp platform.xml $out/
          cp hostfile.txt $out/
        '';
      };
      dc-simgrid-cuda = pkgs-cuda.cudaPackages.backendStdenv.mkDerivation {
        pname = "dc";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = [pkgs.openmpi pkgs.simgrid pkgs-cuda.cudatoolkit];
        buildPhase = "make all BACKEND=simgrid_cuda";
        unpackPhase = ''
          mkdir source
          cp -r --no-preserve=mode,ownership $src/* source/
          cd source
        '';
        installPhase = ''
          mkdir -p $out/bin
          cp bin/dc $out/bin
          cp platform.xml $out/
          cp hostfile.txt $out/
        '';
      };

      dc-cuda = pkgs-cuda.cudaPackages.backendStdenv.mkDerivation {
        pname = "dc-cuda";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = [pkgs.openmpi pkgs-cuda.cudatoolkit akypuera];
        buildPhase = "make all BACKEND=cuda";
        unpackPhase = ''
          mkdir source
          cp -r --no-preserve=mode,ownership $src/* source/
          cd source
        '';
        installPhase = ''
          mkdir -p $out/bin
          cp bin/dc $out/bin/dc
        '';
      };

      rEnv = pkgs.rWrapper.override {
        packages = with pkgs.rPackages; [
          languageserver
          lintr
          here
          digest
          tidyverse
          plotly
        ];
      };
    in {
      devShell = pkgs.mkShell {
        buildInputs = [
          pkgs-cuda.cudatoolkit
          pkgs.openmpi
          pkgs.clang-tools
          pkgs.llvmPackages.openmp
          pkgs.simgrid
          pkgs.vite
          pkgs.pandoc
          rEnv
          akypuera
          pajeng
        ];
        shellHook = ''
          export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          export CUDA_PATH=${pkgs-cuda.cudatoolkit}
        '';
      };
      packages = {
        script = pkgs.writeShellScriptBin "run-dc" ''
          BACKEND=$1

          if [ "$BACKEND" == "cuda" ]; then
            echo "Using CUDA backend..."
            APP_DIR=${dc-cuda}
          elif [ "$BACKEND" == "openmp" ]; then
            echo "Using OpenMP backend..."
            APP_DIR=${dc}
          else
            echo "Error: Missing or invalid argument."
            echo "Usage: nix run .#default -- [openmp|cuda]"
            exit 1
          fi

          size_x=100
          size_y=100
          size_z=100
          absorption=2
          dx=1e-1
          dy=1e-1
          dz=1e-1
          dt=1e-6
          tmax=1e-4

          ${pkgs.openmpi}/bin/mpirun -np 6 --bind-to none $APP_DIR/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc
        '';
        simgrid = pkgs.writeShellScriptBin "run-simgrid" ''
          BACKEND=$1
          SIZE=$2

          if [ "$BACKEND" == "cuda" ]; then
            echo "Using CUDA backend..."
            APP_DIR=${dc-simgrid-cuda}
          elif [ "$BACKEND" == "openmp" ]; then
            echo "Using OpenMP backend..."
            APP_DIR=${dc-simgrid}
          else
            echo "Error: Missing or invalid argument."
            echo "Usage: nix run .#simgrid -- [openmp|cuda]"
            exit 1
          fi

          absorption=2
          size_x=$(( SIZE - 2 * absorption - 8 ))
          size_y=$size_x
          size_z=$size_y
          if (( SIZE % 32 != 0 || SIZE <= 0 )); then
            echo "Error: size must be divisible by 32 and positive"
            exit 1
          fi

          dx=1e-1
          dy=1e-1
          dz=1e-1
          dt=1e-6
          tmax=1e-4

          export OMP_NUM_THREADS=1
          ${pkgs.simgrid}/bin/smpirun -platform $APP_DIR/platform.xml --cfg=smpi/display-timing:yes --cfg=precision/timing:1e-9 --cfg=tracing/precision:9 --cfg=smpi/host-speed:auto -trace --cfg=tracing/filename:dc.trace -hostfile $APP_DIR/hostfile.txt $APP_DIR/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc

          ${pajeng}/bin/pj_dump -l 9 dc.trace | grep ^State > dc.csv
          Rscript ./plot.R
        '';
        comparison = pkgs.writeShellScriptBin "run-dc-comparison" ''
          export PATH=${pkgs.cudatoolkit}/bin:$PATH
          size_x=52
          size_y=52
          size_z=52
          absorption=2
          dx=1e-1
          dy=1e-1
          dz=1e-1
          dt=1e-6
          tmax=1e-4

          echo "Running sequential version to generate ground truth..."
          OMP_NUM_THREADS=1 ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${dc}/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/ground_truth.dc

          echo "Running OpenMP version..."
          ${pkgs.openmpi}/bin/mpirun -np 6 --bind-to none ${dc}/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/openmp_predicted.dc

          echo "Comparing OpenMP with ground truth..."
          mv ./validation/openmp_predicted.dc ./validation/predicted.dc
          Rscript ./validation/CompareResults.R 0

          echo "Running CUDA version..."
          ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${dc-cuda}/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc

          echo "Comparing CUDA with ground truth..."

          echo "--- Comparison with tolerance 1e-3 ---"
          Rscript ./validation/CompareResults.R 1e-3

          echo "--- Comparison with tolerance 1e-5 ---"
          Rscript ./validation/CompareResults.R 1e-5

          echo "--- Comparison with tolerance 1e-7 ---"
          Rscript ./validation/CompareResults.R 1e-7

          echo "--- Comparison with tolerance 0 (exact) ---"
          Rscript ./validation/CompareResults.R 0
        '';
        default = dc;
        cuda = dc-cuda;
      };
      apps = {
        default = {
          type = "app";
          program = "${self.packages.${system}.script}/bin/run-dc";
        };
        simgrid = {
          type = "app";
          program = "${self.packages.${system}.simgrid}/bin/run-simgrid";
        };
        comparison = {
          type = "app";
          program = "${self.packages.${system}.comparison}/bin/run-dc-comparison";
        };
      };
    });
}
