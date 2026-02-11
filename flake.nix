{
  inputs = {
    utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.11";
  };
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
      akypuera = import ./nix/akypuera.nix {inherit pkgs;};
      pajeng = import ./nix/pajeng.nix {inherit pkgs;};
      mpiP = import ./nix/mpiP.nix {inherit pkgs;};

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

      packages = import ./nix/packages.nix {inherit pkgs mpiP akypuera;};
      scripts = import ./nix/scripts.nix {inherit pkgs packages pajeng rEnv akypuera;};
    in {
      devShell = pkgs.mkShell {
        buildInputs =
          [
            pkgs.cudatoolkit
            pkgs.openmpi
            pkgs.clang-tools
            pkgs.llvmPackages.openmp
            pkgs.simgrid
            pkgs.vite
            pkgs.pandoc
            pkgs.pkg-config
            rEnv
            akypuera
            pajeng
          ]
          ++ builtins.attrValues scripts;
        shellHook = ''
          export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          export CUDA_PATH=${pkgs.cudatoolkit}
        '';
      };

      packages =
        packages
        // scripts
        // {
          default = packages.dc-omp-mpip;
        };

      apps = {
        default = {
          type = "app";
          program = "${scripts.run-dc}/bin/run-dc";
        };
        simgrid = {
          type = "app";
          program = "${scripts.run-simgrid}/bin/run-simgrid";
        };
        poti = {
          type = "app";
          program = "${scripts.run-poti-experiments}/bin/run-poti-experiments";
        };
        dc = {
          type = "app";
          program = "${scripts.dc}/bin/dc";
        };
      };
    });
}
