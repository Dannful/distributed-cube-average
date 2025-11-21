{
  inputs = {utils.url = "github:numtide/flake-utils";};
  outputs = {
    self,
    nixpkgs,
    utils,
  }:
    utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {
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
        nativeBuildInputs = with pkgs; [gnumake openmpi mpiP];
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

      dc-cuda = pkgs.stdenv.mkDerivation {
        pname = "dc-cuda";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = with pkgs; [gnumake openmpi mpiP cudatoolkit];
        buildPhase = "make all BACKEND=cuda";
        installPhase = ''
          mkdir -p $out/bin
          cp bin/dc $out/bin/dc
        '';
      };

      mpiP = import ./mpiP.nix {
        inherit pkgs;
      };
      rEnv = pkgs.rWrapper.override {
        packages = with pkgs.rPackages; [
          languageserver
          lintr
          here
          digest
          tidyverse
        ];
      };
    in {
      devShell = pkgs.mkShell {
        buildInputs = [
          pkgs.cudatoolkit
          pkgs.openmpi
          pkgs.clang-tools
          pkgs.llvmPackages.openmp
          pkgs.simgrid
          pkgs.vite
          rEnv
          akypuera
          pajeng
        ];
        shellHook = ''
          export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          export CUDA_PATH=${pkgs.cudatoolkit}
        '';
      };
      packages = {
        script = pkgs.writeShellScriptBin "run-dc" ''
          size_x=100
          size_y=100
          size_z=100
          absorption=2
          dx=1e-1
          dy=1e-1
          dz=1e-1
          dt=1e-6
          tmax=1e-4

          ${pkgs.openmpi}/bin/mpirun -np 8 --bind-to none ${dc}/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc
        '';
        simgrid = pkgs.writeShellScriptBin "run-simgrid" ''
          size_x=52
          size_y=52
          size_z=52
          absorption=2
          dx=1e-1
          dy=1e-1
          dz=1e-1
          dt=1e-6
          tmax=1e-4

          ${pkgs.simgrid}/bin/smpirun -platform ${dc-simgrid}/platform.xml --cfg=smpi/display-timing:yes --cfg=precision/timing:1e-9 --cfg=tracing/precision:9 --cfg=smpi/host-speed:auto -trace --cfg=tracing/filename:dc.trace -hostfile ${dc-simgrid}/hostfile.txt ${dc-simgrid}/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc
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
          ${pkgs.openmpi}/bin/mpirun -np 8 --bind-to none ${dc}/bin/dc --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/openmp_predicted.dc

          echo "Comparing OpenMP with ground truth..."
          mv ./validation/openmp_predicted.dc ./validation/predicted.dc
          Rscript ./validation/CompareResults.R 0

          echo "Running CUDA version..."
          ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${dc-cuda}/bin/dc-cuda --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc

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
