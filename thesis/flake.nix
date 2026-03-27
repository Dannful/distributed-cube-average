{
  inputs = {
    utils.url = "github:numtide/flake-utils";
  };
  outputs = {
    self,
    nixpkgs,
    utils,
  }:
    utils.lib.eachDefaultSystem (
      system: let
        pkgs = nixpkgs.legacyPackages.${system};
        thesis = pkgs.writeShellScriptBin "build-latex" ''
          pdflatex index.tex
          makeglossaries index
          bibtex index.aux
          pdflatex index.tex
          pdflatex index.tex
        '';
      in {
        devShell = pkgs.mkShell {
          buildInputs = with pkgs; [
            texliveFull
          ];
        };
        apps = {
          default = {
            type = "app";
            program = pkgs.lib.getExe thesis;
          };
        };
      }
    );
}
