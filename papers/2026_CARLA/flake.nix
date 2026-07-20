{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.11";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    utils,
  }:
    utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {inherit system;};

      emacsEnv = pkgs.emacs.pkgs.withPackages (epkgs: [
        epkgs.org-contrib # provides ox-extra.el
      ]);

      texliveEnv = pkgs.texlive.combine {
        inherit
          (pkgs.texlive)
          scheme-small
          booktabs
          hyperref
          amsmath
          microtype
          ;
      };

      rEnv = pkgs.rWrapper.override {
        packages = with pkgs.rPackages; [
          ggplot2
          tidyr
          dplyr
        ];
      };

      name = "CARLA_SimGrid_Ghost";

      pdf = pkgs.stdenv.mkDerivation {
        inherit name;
        src = ./.;
        nativeBuildInputs = [emacsEnv texliveEnv pkgs.which];
        buildPhase = ''
          export HOME=$(mktemp -d)
          emacs -batch \
            --eval "(require 'ox-latex)" \
            --eval "(require 'ox-extra)" \
            --eval "(ox-extras-activate '(ignore-headlines))" \
            --eval "(setq org-export-babel-evaluate nil)" \
            --eval "(setq org-confirm-babel-evaluate nil)" \
            --eval "(add-to-list 'org-latex-classes \
              '(\"LLNCS\" \
                \"\\\\documentclass[runningheads]{llncs}\\n [NO-DEFAULT-PACKAGES]\\n [EXTRA]\\n\" \
                (\"\\\\section{%s}\" . \"\\\\section*{%s}\") \
                (\"\\\\subsection{%s}\" . \"\\\\subsection*{%s}\") \
                (\"\\\\subsubsection{%s}\" . \"\\\\subsubsection*{%s}\") \
                (\"\\\\paragraph{%s}\" . \"\\\\paragraph*{%s}\") \
                (\"\\\\subparagraph{%s}\" . \"\\\\subparagraph*{%s}\")))" \
            ${name}.org --funcall org-latex-export-to-latex
          pdflatex ${name}.tex
          bibtex ${name}
          pdflatex ${name}.tex
          pdflatex ${name}.tex
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
        buildInputs = [emacsEnv texliveEnv rEnv];
      };
    });
}
