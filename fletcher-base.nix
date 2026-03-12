{
  pkgs ?
    import <nixpkgs> {
      config = {
        allowUnfree = true;
      };
    },
}:
pkgs.cudaPackages.backendStdenv.mkDerivation {
  pname = "fletcher-base";
  version = "1.0.0";
  src = pkgs.fetchFromGitHub {
    owner = "gabrielfrtg";
    repo = "fletcher-base";
    rev = "55ea0803f5d2afe58be8e6cff8bd9108cae30b43";
    hash = "sha256-heqW2u5mtiQgSif/toVHIKtvPKsRLOokE7n344fx8Ao=";
  };
  nativeBuildInputs = with pkgs; [gcc cudatoolkit];
  buildPhase = ''
    source ./env.sh
    cd original
    make all backend=CUDA
  '';
  installPhase = ''
    mkdir -p $out/bin
    cp ModelagemFletcher.exe $out/bin/fletcher-cuda
  '';
}
