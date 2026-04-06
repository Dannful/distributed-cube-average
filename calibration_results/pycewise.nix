{pkgs ? import <nixpkgs> {}}:
pkgs.python3Packages.buildPythonPackage rec {
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
    maintainers = [];
  };
}
