.PHONY: all clean

# Run make here, or make -C notes from the repository root.
# Requires pandoc and its default LaTeX PDF engine, pdflatex.
# Mermaid fences remain code blocks; diagrams are not rendered automatically.
MARKDOWN_SOURCES := $(wildcard *.md)
GENERATED_PDFS := $(MARKDOWN_SOURCES:.md=.pdf)

all: $(GENERATED_PDFS)

# Rebuild only when the Markdown source or this recipe is newer.
# Map Unicode Delta and right arrow to LaTeX symbols, preserving the source text.
%.pdf: %.md Makefile
	pandoc "$<" \
		-V geometry:margin=1in \
		-V fontsize=11pt \
		-V 'header-includes=\DeclareUnicodeCharacter{0394}{\ensuremath{\Delta}}\DeclareUnicodeCharacter{2192}{\ensuremath{\rightarrow}}' \
		-o "$@"

# Remove only PDFs corresponding to Markdown sources in this directory.
clean:
	rm -f $(GENERATED_PDFS)
