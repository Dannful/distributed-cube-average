{
  description = "SMPI Calibration Environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
  let
    system = "x86_64-linux";
    pkgs = nixpkgs.legacyPackages.${system};
    
    pythonEnv = pkgs.python3.withPackages (ps: with ps; [
      pandas
      numpy
      matplotlib
      jupyter
      scikit-learn
      # pycewise may need to be built from source if not in nixpkgs
    ]);

    rEnv = pkgs.rWrapper.override {
      packages = with pkgs.rPackages; [ Ckmeans_1d_dp ];
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
      ];

      shellHook = ''
        echo "SMPI Calibration Environment Loaded."
        echo "Available tools: smpirun, python, R, cmake."
      '';
    };
  };
}
