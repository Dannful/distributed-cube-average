{
  inputs = { utils.url = "github:numtide/flake-utils"; };
  outputs = { self, nixpkgs, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        distributed-cube-average = pkgs.stdenv.mkDerivation {
          pname = "distributed-cube-average";
          version = "0.1.0";
          src = ./.;
          nativeBuildInputs = with pkgs; [ gnumake openmpi ];
          buildPhase = "make all";
          installPhase = ''
            mkdir -p $out/bin
            cp bin/distributed-cube-average $out/bin/
          '';
        };
      in {
        devShell = pkgs.mkShell {
          buildInputs = with pkgs; [ openmpi clang-tools ];
          shellHook = ''
            export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          '';
        };
        packages.default =
          pkgs.writeShellScriptBin "run-distributed-cube-average" ''
            ${pkgs.openmpi}/bin/mpirun --map-by :OVERSUBSCRIBE -np 27 ${distributed-cube-average}/bin/distributed-cube-average
          '';
        apps.default = {
          type = "app";
          program = "${
              self.packages.${system}.default
            }/bin/run-distributed-cube-average";
        };
      });
}
