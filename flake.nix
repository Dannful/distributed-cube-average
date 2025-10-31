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

      akypuera = import ./akypuera.nix {pkgs = pkgs;};

      pajeng = import ./pajeng.nix {pkgs = pkgs;};

      distributed-cube-average = pkgs.stdenv.mkDerivation {
        pname = "distributed-cube-average";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = with pkgs; [gnumake openmpi akypuera];
        buildPhase = "make all";
        installPhase = ''
          mkdir -p $out/bin
          cp bin/distributed-cube-average $out/bin/
        '';
      };

      distributed-cube-average-cuda = pkgs.stdenv.mkDerivation {
        pname = "distributed-cube-average-cuda";
        version = "0.1.0";
        src = ./.;
        nativeBuildInputs = with pkgs; [gnumake openmpi akypuera cudatoolkit];
        buildPhase = "make all BACKEND=cuda";
        installPhase = ''
          mkdir -p $out/bin
          cp bin/distributed-cube-average $out/bin/
        '';
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
        script = pkgs.writeShellScriptBin "run-distributed-cube-average" ''
          ${pkgs.openmpi}/bin/mpirun -np 8 --bind-to none ${distributed-cube-average}/bin/distributed-cube-average --size-x=100 --size-y=100 --size-z=100 --absorption=6 --dx=2 --dy=3 --dz=4 --dt=0.000110 --time-max=0.00110 --output-file=./validation/predicted.dc
        '';
        comparison = pkgs.writeShellScriptBin "run-distributed-cube-average-comparison" ''
          export PATH=${pkgs.cudatoolkit}/bin:$PATH
          size_x=5
          size_y=5
          size_z=5
          absorption=6
          dx=2
          dy=3
          dz=4
          dt=0.000110
          tmax=0.00110

          echo "Running sequential version to generate ground truth..."
          OMP_NUM_THREADS=1 ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${distributed-cube-average}/bin/distributed-cube-average --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/ground_truth.dc

          echo "Running OpenMP version..."
          ${pkgs.openmpi}/bin/mpirun -np 8 --bind-to none ${distributed-cube-average}/bin/distributed-cube-average --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/openmp_predicted.dc

          echo "Comparing OpenMP with ground truth..."
          mv ./validation/openmp_predicted.dc ./validation/predicted.dc
          Rscript ./validation/CompareResults.R 0

          echo "Running CUDA version..."
          ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${distributed-cube-average-cuda}/bin/distributed-cube-average --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/cuda_predicted.dc

          echo "Comparing CUDA with ground truth..."
          mv ./validation/cuda_predicted.dc ./validation/predicted.dc

          echo "--- Comparison with tolerance 1e-3 ---"
          Rscript ./validation/CompareResults.R 1e-3

          echo "--- Comparison with tolerance 1e-5 ---"
          Rscript ./validation/CompareResults.R 1e-5

          echo "--- Comparison with tolerance 1e-7 ---"
          Rscript ./validation/CompareResults.R 1e-7

          echo "--- Comparison with tolerance 0 (exact) ---"
          Rscript ./validation/CompareResults.R 0
        '';
        default = distributed-cube-average;
        cuda = distributed-cube-average-cuda;
      };
      apps.default = {
        type = "app";
        program = "${self.packages.${system}.script}/bin/run-distributed-cube-average";
      };
      apps.cuda = {
        type = "app";
        program = "${self.packages.${system}.script}/bin/run-distributed-cube-average";
      };
      apps.comparison = {
        type = "app";
        program = "${self.packages.${system}.comparison}/bin/run-distributed-cube-average-comparison";
      };
    });
}
