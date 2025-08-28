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
        rEnv = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [ languageserver lintr here digest ];
        };

        akypuera = pkgs.stdenv.mkDerivation rec {
          pname = "akypuera";
          version = "0.1.0";

          src = pkgs.fetchFromGitHub {
            owner = "schnorr";
            repo = "akypuera";
            rev = "c78b02ec4d154a1a60730bcce7e6c4cd23f43604";
            sha256 = "sha256-pljg0tlIj4RKN/DGUTeeZu9+rvBrieIfz2HwRmXfkUI=";
            fetchSubmodules = true;
          };

          nativeBuildInputs = with pkgs; [ cmake gnumake gcc autoPatchelfHook ];

          buildInputs = with pkgs; [ openmpi ];

          buildPhase = ''
            runHook preBuild

            mkdir -p build
            cd build
            cmake ..
            make

            runHook postBuild
          '';

          installPhase = ''
            runHook preInstall

            cp libaky.so $out/

            runHook postInstall
          '';

          meta = with pkgs.lib; {
            description = "A Nix derivation for ${pname}";
            homepage = "https://github.com/${src.owner}/${src.repo}";
            license = licenses.gpl3;
          };
        };
      in {
        devShell = pkgs.mkShell {
          buildInputs = [ pkgs.openmpi pkgs.clang-tools rEnv ];
          shellHook = ''
            export PATH=${pkgs.clang-tools}/bin/clangd:$PATH
          '';
        };
        packages.default =
          pkgs.writeShellScriptBin "run-distributed-cube-average" ''
            LD_PRELOAD=${akypuera}/libaky.so:${pkgs.openmpi}/lib/libmpi.so.40 ${pkgs.openmpi}/bin/mpirun --map-by :OVERSUBSCRIBE -np 8 ${distributed-cube-average}/bin/distributed-cube-average 4 4 4 2 1 ./predicted.dc
          '';
        apps.default = {
          type = "app";
          program = "${
              self.packages.${system}.default
            }/bin/run-distributed-cube-average";
        };
      });
}
