TEX_DIR    := tex
BUILD_DIR  := build
PDF_DIR    := $(BUILD_DIR)/pdf
DOCX_DIR   := $(BUILD_DIR)/docx

SOURCES    := $(wildcard $(TEX_DIR)/*.tex)
PDFS       := $(patsubst $(TEX_DIR)/%.tex,$(PDF_DIR)/%.pdf,$(SOURCES))
DOCXS      := $(patsubst $(TEX_DIR)/%.tex,$(DOCX_DIR)/%.docx,$(SOURCES))

# Three-paper split (tex/papers/<part>/<paper>.tex); built in place with tectonic.
PAPERS_DIR    := $(TEX_DIR)/papers
PAPER_SOURCES := $(wildcard $(PAPERS_DIR)/*/*.tex)
PAPER_PDFS    := $(PAPER_SOURCES:.tex=.pdf)

LATEXMK    ?= latexmk
LATEX_ENG  ?= lualatex
PANDOC     ?= pandoc
TECTONIC   ?= tectonic

LATEXMK_FLAGS := -$(LATEX_ENG) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=$(PDF_DIR)

.PHONY: all pdf docx papers clean distclean

all: pdf docx papers

pdf: $(PDFS)

docx: $(DOCXS)

papers: $(PAPER_PDFS)

$(PAPERS_DIR)/%.pdf: $(PAPERS_DIR)/%.tex
	cd $(dir $<) && $(TECTONIC) -X compile $(notdir $<)

$(PDF_DIR)/%.pdf: $(TEX_DIR)/%.tex | $(PDF_DIR)
	$(LATEXMK) $(LATEXMK_FLAGS) $<

$(DOCX_DIR)/%.docx: $(TEX_DIR)/%.tex | $(DOCX_DIR)
	$(PANDOC) $< -o $@

$(PDF_DIR) $(DOCX_DIR):
	mkdir -p $@

clean:
	$(LATEXMK) -C -output-directory=$(PDF_DIR) $(SOURCES) 2>/dev/null || true
	rm -f $(PDF_DIR)/*.aux $(PDF_DIR)/*.log $(PDF_DIR)/*.out $(PDF_DIR)/*.fls $(PDF_DIR)/*.fdb_latexmk

distclean:
	rm -rf $(BUILD_DIR)
	rm -f $(PAPER_PDFS)
