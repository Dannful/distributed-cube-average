# Makefile.nix — compila o PDF via Nix Flake
#
# Uso:
#   make -f Makefile.nix        # compila o PDF (nix build)
#   make -f Makefile.nix shell  # entra no devShell (emacs + texlive + R)
#   make -f Makefile.nix clean  # limpa o store do projeto (nix flake cleanup opcional)

NAME := CARLA_SimGrid_Ghost

.PHONY: all pdf shell clean

all: pdf

# Compila o PDF usando o flake.nix
pdf:
	nix build .#pdf
	@echo "PDF gerado: $(NAME).pdf"
	@echo "Localizado em: $(shell nix eval --raw .#pdf)/$(NAME).pdf"

# Abre um shell de desenvolvimento com emacs, texlive e R
shell:
	nix develop

# Limpeza — remove o lock e reconstrói
clean:
	nix flake lock --update-input nixpkgs
	@echo "flake.lock atualizado. Execute 'make -f Makefile.nix pdf' para reconstruir."
