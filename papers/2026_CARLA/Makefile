NAME = CARLA_SimGrid_Ghost

EMACS_SETUP = \
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
      (\"\\\\subparagraph{%s}\" . \"\\\\subparagraph*{%s}\")))"

.PHONY: all clean distclean

all: $(NAME).pdf

$(NAME).tex: $(NAME).org
	emacs -batch $(EMACS_SETUP) $< --funcall org-latex-export-to-latex

$(NAME).pdf: $(NAME).tex biblio.bib $(wildcard img/*)
	pdflatex $(NAME).tex
	bibtex $(NAME)
	pdflatex $(NAME).tex
	pdflatex $(NAME).tex

clean:
	rm -f $(NAME).aux $(NAME).bbl $(NAME).blg $(NAME).log $(NAME).out $(NAME).tex *~

distclean: clean
	rm -f $(NAME).pdf
