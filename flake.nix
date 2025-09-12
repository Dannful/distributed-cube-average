{
  inputs = { utils.url = "github:numtide/flake-utils"; };
  outputs = { self, nixpkgs, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

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

          cmakeFlags = [ "DC_MAKE_SKIP_BUILD_RPATH=ON" ];

          nativeBuildInputs = with pkgs; [ cmake gnumake gcc ];

          buildInputs = with pkgs; [ openmpi gcc ];

          installPhase = ''
            runHook preInstall

            mkdir -p $out/bin
            mkdir -p $out/lib

            cp aky_converter $out/bin/
            cp libaky.so $out/lib/
            cp librastro/librastro.so $out/lib/
            cp poti/libpoti.so.8 $out/lib

            local RPATH="${pkgs.lib.makeLibraryPath buildInputs}:$out/lib"

            patchelf --set-rpath "$RPATH" $out/lib/libaky.so
            patchelf --set-rpath "$RPATH" $out/bin/aky_converter

            runHook postInstall
          '';

          meta = with pkgs.lib; {
            description = "A Nix derivation for ${pname}";
            homepage = "https://github.com/${src.owner}/${src.repo}";
            license = licenses.gpl3;
          };
        };

        pajeng = pkgs.stdenv.mkDerivation rec {
          pname = "pajeng";
          version = "1.3.10";

          src = pkgs.fetchFromGitHub {
            owner = "schnorr";
            repo = "pajeng";
            rev = "87d2d263020339defddeeaf6586723d0f840f5a4";
            sha256 = "sha256-g3aT5SNwrhk7d/f5ElJfSqXQ6MFnsuVQ0fSpbiH94Y0=";
          };

          nativeBuildInputs = with pkgs; [
            cmake
            gnumake
            gcc
            flex
            bison
            asciidoc
            boost
            fmt
          ];

          installPhase = ''
            runHook preInstall

            mkdir -p $out/bin
            mkdir -p $out/lib

            cp pj_dump $out/bin/
            cp pj_equals $out/bin/
            cp src/libpaje/libpaje.so.2 $out/lib

            local RPATH="${pkgs.lib.makeLibraryPath buildInputs}:$out/lib"

            patchelf --set-rpath "$RPATH" $out/bin/pj_dump
            patchelf --set-rpath "$RPATH" $out/bin/pj_equals

            runHook postInstall
          '';

          buildInputs = with pkgs; [ openmpi gcc ];

          meta = with pkgs.lib; {
            description = "A Nix derivation for ${pname}";
            homepage = "https://github.com/${src.owner}/${src.repo}";
            license = licenses.gpl3;
          };
        };

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
            ${pkgs.openmpi}/bin/mpirun --map-by :OVERSUBSCRIBE -np 2 ${distributed-cube-average}/bin/distributed-cube-average --size-x=4 --size-y=4 --size-z=4 --dx=1 --dy=1 --dz=1 --dt=1 --time-max=6 --stencil-size=1 --absorption=3 --output-file=./validation/predicted.dc
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
