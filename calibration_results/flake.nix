{
  description = "SMPI Calibration Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = {
    self,
    nixpkgs,
  }: let
    system = "x86_64-linux";
    pkgs = nixpkgs.legacyPackages.${system};

    pycewise = pkgs.python3Packages.buildPythonPackage rec {
      pname = "pycewise";
      version = "unstable-2021-05-18";

      src = pkgs.fetchFromGitHub {
        owner = "Ezibenroc";
        repo = "pycewise";
        rev = "7f7608e";
        hash = "sha256-FTg9GTLC0if4bWiFTb0DLY6bcgYa1R5bGcKBOVZ3nK0=";
      };

      format = "setuptools";

      postPatch = ''
        substituteInPlace setup.py \
          --replace-fail "'__git_version__': git_version()," "'__git_version__': '${version}',"
      '';

      propagatedBuildInputs = with pkgs.python3Packages; [
        numpy
        statsmodels
        matplotlib
        graphviz
        palettable
      ];

      pythonImportsCheck = ["pycewise"];

      doCheck = false;

      meta = with pkgs.lib; {
        description = "Python module to compute a segmented linear regression";
        homepage = "https://github.com/Ezibenroc/pycewise";
        license = licenses.mit;
      };
    };

    platform-calibration = pkgs.stdenv.mkDerivation {
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
    };

    pythonEnv = pkgs.python3.withPackages (ps:
      with ps; [
        pandas
        numpy
        matplotlib
        jupyter
        scikit-learn
        plotnine
        pycewise
        papermill
        pyyaml
      ]);

    rEnv = pkgs.rWrapper.override {
      packages = with pkgs.rPackages; [Ckmeans_1d_dp];
    };
  in {
    devShells.${system}.default = pkgs.mkShell {
      buildInputs = [
        pkgs.simgrid
        pkgs.openmpi
        pkgs.cmake
        pkgs.gcc
        pythonEnv
        rEnv
        platform-calibration
      ];

      shellHook = ''
        echo "SMPI Calibration Environment Loaded."
        echo "Available tools: smpirun, python, R, cmake, calibrate, bp_search1, bp_search2"
      '';
    };
  };
}
