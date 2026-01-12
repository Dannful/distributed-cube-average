{ pkgs ? import <nixpkgs> { } }:

pkgs.stdenv.mkDerivation rec {
  pname = "pajeng";
  version = "1.3.10";

  src = pkgs.fetchFromGitHub {
    owner = "schnorr";
    repo = "pajeng";
    rev = "87d2d263020339defddeeaf6586723d0f840f5a4";
    sha256 = "sha256-g3aT5SNwrhk7d/f5ElJfSqXQ6MFnsuVQ0fSpbiH94Y0=";
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
    cp pj_equals $out/bin/
    cp src/libpaje/libpaje.so.2 $out/lib

    local gcc_lib=$(gcc -print-file-name=libstdc++.so)
    local RPATH="${
      pkgs.lib.makeLibraryPath buildInputs
    }:$out/lib:$(dirname $gcc_lib)"

    # Set the RPATH on the binaries
    patchelf --set-rpath "$RPATH" $out/bin/pj_dump
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
