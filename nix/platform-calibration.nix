{pkgs ? import <nixpkgs> {}}:
pkgs.stdenv.mkDerivation {
  pname = "platform-calibration";
  version = "1.0.0";

  src = pkgs.fetchgit {
    url = "https://framagit.org/simgrid/platform-calibration.git";
    rev = "ec784f1f19e51e2d759ba80553efd1c2fed34148";
    hash = "sha256-24cmhrpoSm1iKo9gRdf372iBOIiHyPPluUAh3gyvHWQ=";
  };

  nativeBuildInputs = with pkgs; [openmpi];

  buildPhase = ''
    cd src/calibration
    make
  '';

  installPhase = ''
    mkdir -p $out/bin
    cp calibrate $out/bin
    cp bp_search1 $out/bin
    cp bp_search2 $out/bin
  '';

  buildInputs = with pkgs; [openmpi];

  meta = with pkgs.lib; {
    description = "Procedures and results of calibration efforts for SimGrid platforms";
    homepage = "https://framagit.org/simgrid/platform-calibration";
    license = licenses.mit;
    platforms = platforms.unix;
    maintainers = [];
  };
}
