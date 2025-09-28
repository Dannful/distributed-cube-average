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
          packages = with pkgs.rPackages; [ languageserver lintr here digest ];
        };

      in {
        devShell = pkgs.mkShell {
          buildInputs = [ pkgs.openmpi pkgs.clang-tools rEnv akypuera pajeng ];
          shellHook = ''
            export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          '';
        };
        packages.script =
          pkgs.writeShellScriptBin "run-distributed-cube-average" ''
            ${pkgs.openmpi}/bin/mpirun -np 8 ${distributed-cube-average}/bin/distributed-cube-average --size-x=20 --size-y=20 --size-z=20 --absorption=6 --dx=2 --dy=3 --dz=4 --dt=0.000110 --time-max=3 --output-file=./validation/predicted.dc
          '';
        packages.default = distributed-cube-average;
        apps.default = {
          type = "app";
          program = "${
              self.packages.${system}.script
            }/bin/run-distributed-cube-average";
        };
      });
}
