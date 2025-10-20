{
  inputs = { utils.url = "github:numtide/flake-utils"; };
  outputs = { self, nixpkgs, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        akypuera = import ./akypuera.nix { pkgs = pkgs; };

        pajeng = import ./pajeng.nix { pkgs = pkgs; };

        distributed-cube-average = pkgs.stdenv.mkDerivation {
          pname = "distributed-cube-average";
          version = "0.1.0";
          src = ./.;
          nativeBuildInputs = with pkgs; [ gnumake openmpi akypuera ];
          buildPhase = "make all";
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
            pkgs.openmpi
            pkgs.clang-tools
            pkgs.llvmPackages.openmp
            rEnv
            akypuera
            pajeng
          ];
          shellHook = ''
            export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          '';
        };
        packages.script =
          pkgs.writeShellScriptBin "run-distributed-cube-average" ''
            ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${distributed-cube-average}/bin/distributed-cube-average --size-x=100 --size-y=100 --size-z=100 --absorption=6 --dx=2 --dy=3 --dz=4 --dt=0.000110 --time-max=1 --output-file=./validation/predicted.dc
          '';
        packages.comparison =
          pkgs.writeShellScriptBin "run-distributed-cube-average-comparison" ''
            size_x=5
            size_y=5
            size_z=5
            absorption=6
            dx=2
            dy=3
            dz=4
            dt=0.000110
            tmax=0.00110
            OMP_NUM_THREADS=1 ${pkgs.openmpi}/bin/mpirun -np 1 --bind-to none ${distributed-cube-average}/bin/distributed-cube-average --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/ground_truth.dc
            ${pkgs.openmpi}/bin/mpirun -np 8 --bind-to none ${distributed-cube-average}/bin/distributed-cube-average --size-x=$size_x --size-y=$size_y --size-z=$size_z --absorption=$absorption --dx=$dx --dy=$dy --dz=$dz --dt=$dt --time-max=$tmax --output-file=./validation/predicted.dc
            Rscript ./validation/CompareResults.R
          '';
        packages.default = distributed-cube-average;
        apps.default = {
          type = "app";
          program = "${
              self.packages.${system}.script
            }/bin/run-distributed-cube-average";
        };
        apps.comparison = {
          type = "app";
          program = "${
              self.packages.${system}.comparison
            }/bin/run-distributed-cube-average-comparison";
        };
      });
}
