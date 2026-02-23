# Minimal makefile for Sphinx documentation

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help html latexpdf latexpdf-ja latexpdf-en clean

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)

latexpdf:
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	$(MAKE) -C "$(BUILDDIR)/latex" all-pdf

latexpdf-ja:
	SPHINX_LANG=ja $(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex-ja" $(SPHINXOPTS) $(O)
	$(MAKE) -C "$(BUILDDIR)/latex-ja" all-pdf

latexpdf-en:
	SPHINX_LANG=en $(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex-en" $(SPHINXOPTS) $(O)
	$(MAKE) -C "$(BUILDDIR)/latex-en" all-pdf

clean:
	rm -rf $(BUILDDIR)
