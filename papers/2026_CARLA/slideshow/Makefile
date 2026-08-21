NAME = slides

EMACS_SETUP = \
  --eval "(require 'ox-beamer)" \
  --eval "(require 'ox-extra)" \
  --eval "(ox-extras-activate '(ignore-headlines))" \
  --eval "(setq org-export-babel-evaluate nil)" \
  --eval "(setq org-confirm-babel-evaluate nil)"

.PHONY: all clean distclean

all: $(NAME).pdf

$(NAME).tex: $(NAME).org
	emacs -batch $(EMACS_SETUP) $< --funcall org-beamer-export-to-latex

# Figures are referenced relative to the parent paper directory (../img).
$(NAME).pdf: $(NAME).tex $(wildcard ../img/*)
	pdflatex -interaction=nonstopmode -halt-on-error $(NAME).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(NAME).tex

clean:
	rm -f $(NAME).aux $(NAME).log $(NAME).nav $(NAME).out $(NAME).snm \
	      $(NAME).toc $(NAME).tex $(NAME)*.vrb *~

distclean: clean
	rm -f $(NAME).pdf
