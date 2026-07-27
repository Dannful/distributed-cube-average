NAME = SBAC_Fletcher_MPI
ARCHIVE = replicable
ARCHIVE_FILES = $(NAME).org refs.bib Makefile img/*

EMACS_SETUP = \
  --eval "(require 'ox-latex)" \
  --eval "(require 'ox-extra)" \
  --eval "(ox-extras-activate '(ignore-headlines))" \
  --eval "(setq org-export-babel-evaluate nil)" \
  --eval "(setq org-confirm-babel-evaluate nil)" \
  --eval "(add-to-list 'org-latex-classes \
    '(\"IEEEtran\" \
      \"\\\\documentclass{IEEEtran}\" \
      (\"\\\\section{%s}\" . \"\\\\section*{%s}\") \
      (\"\\\\subsection{%s}\" . \"\\\\subsection*{%s}\") \
      (\"\\\\subsubsection{%s}\" . \"\\\\subsubsection*{%s}\") \
      (\"\\\\paragraph{%s}\" . \"\\\\paragraph*{%s}\") \
      (\"\\\\subparagraph{%s}\" . \"\\\\subparagraph*{%s}\")))"

.PHONY: all clean distclean

all: $(NAME).pdf

refs.bib: $(NAME).org
	emacs -batch \
	  --eval "(setq org-confirm-babel-evaluate nil)" \
	  $< --funcall org-babel-tangle

$(NAME).tex: $(NAME).org refs.bib
	emacs -batch $(EMACS_SETUP) $< --funcall org-latex-export-to-latex

$(NAME).pdf: $(NAME).tex refs.bib $(wildcard img/*)
	pdflatex $(NAME).tex
	bibtex $(NAME)
	pdflatex $(NAME).tex
	pdflatex $(NAME).tex

$(ARCHIVE).tgz: $(ARCHIVE_FILES)
	tar --xform "s|^|$(ARCHIVE)/|" -zcf $@ $^

clean:
	rm -f $(NAME).aux $(NAME).bbl $(NAME).blg $(NAME).log $(NAME).out $(NAME).tex refs.bib *~

distclean: clean
	rm -f $(NAME).pdf
