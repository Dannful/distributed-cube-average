{pkgs ? import <nixpkgs>}:
pkgs.stdenv.mkDerivation {
  pname = "mpiP";

  version = "3.5";

  src = pkgs.fetchFromGitHub {
    owner = "LLNL";

    repo = "mpiP";

    rev = "8ff38c37777111543307fa40274caa96be8a916b";

    sha256 = "sha256-BTp+Kx3GXwqlN1TBBXxcyvMQmdmC0M6rrnEdJKfeIPE=";
  };

  nativeBuildInputs = [pkgs.python3 pkgs.openmpi pkgs.llvmPackages.libunwind];
}
