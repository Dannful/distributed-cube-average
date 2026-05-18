{ pkgs ? import <nixpkgs> { } }:

pkgs.stdenv.mkDerivation rec {
  pname = "pajeng";
  version = "1.3.10";

  src = pkgs.fetchFromGitHub {
    owner = "schnorr";
    repo = "pajeng";
    rev = "e8d14fcc7a5fbec28d7458e8b856a45c545cb9be";
    sha256 = "sha256-s3us1DbONxWhSYuSe6Jay23q0ZgI2ZegE3SKgoAIEWA=";
  };

  nativeBuildInputs = [
    pkgs.cmake
    pkgs.gnumake
    pkgs.gcc
    pkgs.flex
    pkgs.bison
    pkgs.asciidoc
    pkgs.boost
    pkgs.fmt
  ];

  installPhase = ''
    runHook preInstall

    mkdir -p $out/bin
    mkdir -p $out/lib

    cp pj_dump $out/bin/
    cp pj_dump_csv $out/bin/
    cp pj_equals $out/bin/
    cp src/libpaje/libpaje.so.2 $out/lib

    local gcc_lib=$(gcc -print-file-name=libstdc++.so)
    local RPATH="${
      pkgs.lib.makeLibraryPath buildInputs
    }:$out/lib:$(dirname $gcc_lib)"

    # Set the RPATH on the binaries
    patchelf --set-rpath "$RPATH" $out/bin/pj_dump
    patchelf --set-rpath "$RPATH" $out/bin/pj_dump_csv
    patchelf --set-rpath "$RPATH" $out/bin/pj_equals

    runHook postInstall
  '';

  buildInputs = [ pkgs.gcc pkgs.flex pkgs.fmt ];

  meta = with pkgs.lib; {
    description = "A Nix derivation for ${pname}";
    homepage = "https://github.com/${src.owner}/${src.repo}";
    license = licenses.gpl3;
  };
}
