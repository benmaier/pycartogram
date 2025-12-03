# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
NOTEBOOKDIR   = ../notebooks

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile notebooks clean-notebooks

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Copy notebooks from ../notebooks to examples/
notebooks:
	@echo "Copying notebooks from $(NOTEBOOKDIR) to examples/..."
	@mkdir -p examples
	@if [ -d "$(NOTEBOOKDIR)" ] && [ "$$(ls -A $(NOTEBOOKDIR)/*.ipynb 2>/dev/null)" ]; then \
		cp $(NOTEBOOKDIR)/*.ipynb examples/; \
		echo "Copied $$(ls $(NOTEBOOKDIR)/*.ipynb | wc -l | tr -d ' ') notebook(s)"; \
	else \
		echo "No notebooks found in $(NOTEBOOKDIR)"; \
	fi

# Clean copied notebooks
clean-notebooks:
	@echo "Removing copied notebooks from examples/..."
	@rm -f examples/*.ipynb

# Build HTML with notebooks
html: notebooks
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Custom targets
livehtml: notebooks
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)

clean: clean-notebooks
	rm -rf $(BUILDDIR)/*
