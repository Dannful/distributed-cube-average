{
  pkgs,
  mpiP,
  akypuera,
}: let
  mkDc = {
    backend,
    profile,
    name,
    extraBuildInputs ? [],
    extraNativeBuildInputs ? [],
  }:
    pkgs.stdenv.mkDerivation {
      pname = name;
      version = "0.1.0";
      src = ../.;
      nativeBuildInputs = with pkgs; [gnumake] ++ extraNativeBuildInputs;
      buildInputs = with pkgs; [openmpi] ++ extraBuildInputs;
      buildPhase = "make all BACKEND=${backend} PROFILE=${profile}";
      installPhase = ''
        mkdir -p $out/bin
        cp bin/dc $out/bin
      '';
    };

  mkDcCuda = {
    backend,
    profile,
    name,
    extraBuildInputs ? [],
    extraNativeBuildInputs ? [],
  }:
    pkgs.cudaPackages.backendStdenv.mkDerivation {
      pname = name;
      version = "0.1.0";
      src = ../.;
      nativeBuildInputs = [pkgs.cudatoolkit] ++ extraNativeBuildInputs;
      buildInputs = [pkgs.openmpi] ++ extraBuildInputs;
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
in {
  # OpenMP variants
  dc-omp-mpip = mkDc {
    name = "dc-omp-mpip";
    backend = "openmp";
    profile = "mpip";
    extraBuildInputs = [mpiP];
  };

  dc-omp-aky = mkDc {
    name = "dc-omp-aky";
    backend = "openmp";
    profile = "akypuera";
    extraBuildInputs = [akypuera];
  };

  # CUDA variants
  dc-cuda-mpip = mkDcCuda {
    name = "dc-cuda-mpip";
    backend = "cuda";
    profile = "mpip";
    extraBuildInputs = [mpiP];
  };

  dc-cuda-aky = mkDcCuda {
    name = "dc-cuda-aky";
    backend = "cuda";
    profile = "akypuera";
    extraBuildInputs = [akypuera];
  };

  # SimGrid variants
  dc-simgrid = pkgs.stdenv.mkDerivation {
    pname = "dc-simgrid";
    version = "0.1.0";
    src = ../.;
    nativeBuildInputs = with pkgs; [gnumake simgrid];
    buildInputs = with pkgs; [openmpi];
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
    nativeBuildInputs = with pkgs; [simgrid cudatoolkit];
    buildInputs = with pkgs; [openmpi];
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

  dc-simgrid-platform = pkgs.stdenv.mkDerivation {
    pname = "dc-simgrid-platform";
    version = "0.1.0";
    src = ../simgrid-config;
    nativeBuildInputs = with pkgs; [ pkg-config ];
    buildInputs = with pkgs; [ simgrid ];
    
    buildPhase = ''
      SG_CFLAGS=$(pkg-config --cflags simgrid)
      SG_LIBS=$(pkg-config --libs simgrid)
      
      echo "Compiling Shared Object (libplatform.so)..."
      g++ -std=c++17 -shared -fPIC -o libplatform.so platform_s4u.cpp $SG_CFLAGS $SG_LIBS
      
      echo "Compiling Generator Executable (generate_artifacts)..."
      g++ -std=c++17 -o generate_artifacts platform_s4u.cpp $SG_CFLAGS $SG_LIBS
    '';
    
    installPhase = ''
      mkdir -p $out/lib $out/bin
      cp libplatform.so $out/lib/
      cp generate_artifacts $out/bin/
    '';
  };
}
