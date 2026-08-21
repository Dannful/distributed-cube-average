{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.11";
    utils.url = "github:numtide/flake-utils";
    # Parent paper flake, used only to obtain the figures under img/.
    paper.url = "path:..";
  };

  outputs = {
    self,
    nixpkgs,
    utils,
    paper,
  }:
    utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {inherit system;};

      emacsEnv = pkgs.emacs.pkgs.withPackages (epkgs: [
        epkgs.org-contrib
      ]);

      texliveEnv = pkgs.texlive.withPackages (ps:
        with ps; [
          scheme-small
          collection-latexrecommended
          collection-latexextra
          collection-fontsrecommended
          beamer
          ifsym # required by beamerthemeInf.sty
        ]);

      name = "slides";

      pdf = pkgs.stdenv.mkDerivation {
        inherit name;
        src = ./.;
        nativeBuildInputs = [emacsEnv texliveEnv pkgs.which];
        buildPhase = ''
          export HOME=$(mktemp -d)
          # Figures are referenced relative to the parent (../img);
          # stage them one level above the source root.
          cp -rL ${paper}/img ../img
          emacs -batch \
            --eval "(require 'ox-beamer)" \
            --eval "(require 'ox-extra)" \
            --eval "(ox-extras-activate '(ignore-headlines))" \
            --eval "(setq org-export-babel-evaluate nil)" \
            --eval "(setq org-confirm-babel-evaluate nil)" \
            ${name}.org --funcall org-beamer-export-to-latex
          pdflatex -interaction=nonstopmode -halt-on-error ${name}.tex
          pdflatex -interaction=nonstopmode -halt-on-error ${name}.tex
        '';
        installPhase = ''
          install -Dm644 ${name}.pdf $out/${name}.pdf
        '';
      };
    in {
      packages = {
        default = pdf;
        inherit pdf;
      };

      devShells.default = pkgs.mkShell {
        buildInputs = [emacsEnv texliveEnv];
      };
    });
}
