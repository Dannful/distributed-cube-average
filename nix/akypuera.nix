{pkgs ? import <nixpkgs> {}}:
pkgs.stdenv.mkDerivation rec {
  pname = "akypuera";
  version = "0.1.0";

  src = pkgs.fetchFromGitHub {
    owner = "schnorr";
    repo = "akypuera";
    rev = "c78b02ec4d154a1a60730bcce7e6c4cd23f43604";
    sha256 = "sha256-pljg0tlIj4RKN/DGUTeeZu9+rvBrieIfz2HwRmXfkUI=";
    fetchSubmodules = true;
  };

  nativeBuildInputs = [pkgs.cmake pkgs.gnumake pkgs.gcc];

  cmakeFlags = ["-DCMAKE_POLICY_VERSION_MINIMUM=3.5"];

  buildInputs = [pkgs.openmpi pkgs.gcc];

  installPhase = ''
    runHook preInstall

    mkdir -p $out/bin
    mkdir -p $out/lib

    cp aky_converter $out/bin/
    cp libaky.so $out/lib/
    cp librastro/librastro.so $out/lib/
    cp poti/libpoti.so.8 $out/lib

    local gcc_lib=$(gcc -print-file-name=libstdc++.so)
    local RPATH="${
      pkgs.lib.makeLibraryPath buildInputs
    }:$out/lib:$(dirname $gcc_lib)"

    patchelf --set-rpath "$RPATH" $out/lib/libaky.so
    patchelf --set-rpath "$RPATH" $out/bin/aky_converter

    runHook postInstall
  '';

  meta = with pkgs.lib; {
    description = "A Nix derivation for ${pname}";
    homepage = "https://github.com/${src.owner}/${src.repo}";
    license = licenses.gpl3;
  };
}
