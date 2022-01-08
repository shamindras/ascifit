# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME = `sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
PKGVERS = `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`

.PHONY : conda_cb01 conda_rem_cb01 clean

all: check

build: install_deps
	R CMD build .

check: build
	R CMD check --no-manual $(PKGNAME)_$(PKGVERS).tar.gz

install_deps:
	Rscript \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE)'

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

conda_cb01:
	conda update -n base -c defaults conda
	conda env create -f=./conda_envs/cb01/environment.yml

conda_rem_cb01:
	conda activate base
	conda remove --name cb01 --all

clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck
	rm -f reading_list.pdf *.out *.aux
	rm -f *.bbl *.blg *.log *.toc *.ptb *.tod *.fls *.fdb_latexmk *.lof *.rds .DS_Store
	find . -name '.DS_Store' -type f -delete
	rm -rf "#README.md#"
	rm -rf ".#README.md"
	rm -rf "#Makefile#"
	rm -rf "#.Makefile"
