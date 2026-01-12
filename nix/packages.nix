{ pkgs, mpiP, akypuera }:
let
  mkDc = { backend, profile, name, extraNativeBuildInputs ? [], extraLdFlags ? "" }:
    pkgs.stdenv.mkDerivation {
      pname = name;
      version = "0.1.0";
      src = ../.;
      nativeBuildInputs = with pkgs; [gnumake openmpi] ++ extraNativeBuildInputs;
      buildPhase = "make all BACKEND=${backend} PROFILE=${profile}";
      installPhase = ''
        mkdir -p $out/bin
        cp bin/dc $out/bin
      '';
    };
  
  mkDcCuda = { backend, profile, name, extraNativeBuildInputs ? [], extraLdFlags ? "" }:
    pkgs.cudaPackages.backendStdenv.mkDerivation {
      pname = name;
      version = "0.1.0";
      src = ../.;
      nativeBuildInputs = [pkgs.openmpi pkgs.cudatoolkit] ++ extraNativeBuildInputs;
      buildPhase = "make all BACKEND=${backend} PROFILE=${profile}";
      unpackPhase = ''
        mkdir source
        cp -r --no-preserve=mode,ownership $src/* source/
        cd source
      '';
      installPhase = ''
        mkdir -p $out/bin
        cp bin/dc $out/bin/dc
      '';
    };
in
{
  # OpenMP variants
  dc-omp-mpip = mkDc {
    name = "dc-omp-mpip";
    backend = "openmp";
    profile = "mpip";
    extraNativeBuildInputs = [ mpiP ];
  };

  dc-omp-aky = mkDc {
    name = "dc-omp-aky";
    backend = "openmp";
    profile = "akypuera";
    extraNativeBuildInputs = [ akypuera ];
  };

  # CUDA variants
  dc-cuda-mpip = mkDcCuda {
    name = "dc-cuda-mpip";
    backend = "cuda";
    profile = "mpip";
    extraNativeBuildInputs = [ mpiP ];
  };

  dc-cuda-aky = mkDcCuda {
    name = "dc-cuda-aky";
    backend = "cuda";
    profile = "akypuera";
    extraNativeBuildInputs = [ akypuera ];
  };

  # SimGrid variants (unchanged)
  dc-simgrid = pkgs.stdenv.mkDerivation {
    pname = "dc-simgrid";
    version = "0.1.0";
    src = ../.;
    nativeBuildInputs = with pkgs; [gnumake openmpi simgrid];
    buildPhase = "make all BACKEND=simgrid";
    installPhase = ''
      mkdir -p $out/bin
      cp bin/dc $out/bin
      cp platform.xml $out/
      cp hostfile.txt $out/
    '';
  };

  dc-simgrid-cuda = pkgs.cudaPackages.backendStdenv.mkDerivation {
    pname = "dc-simgrid-cuda";
    version = "0.1.0";
    src = ../.;
    nativeBuildInputs = with pkgs; [openmpi simgrid cudatoolkit];
    buildPhase = "make all BACKEND=simgrid_cuda";
    unpackPhase = ''
      mkdir source
      cp -r --no-preserve=mode,ownership $src/* source/
      cd source
    '';
    installPhase = ''
      mkdir -p $out/bin
      cp bin/dc $out/bin
      cp platform.xml $out/
      cp hostfile.txt $out/
    '';
  };
}