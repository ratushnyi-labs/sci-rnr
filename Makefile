TEX_DIR    := tex
BUILD_DIR  := build
PDF_DIR    := $(BUILD_DIR)/pdf
DOCX_DIR   := $(BUILD_DIR)/docx

SOURCES    := $(wildcard $(TEX_DIR)/*.tex)
PDFS       := $(patsubst $(TEX_DIR)/%.tex,$(PDF_DIR)/%.pdf,$(SOURCES))
DOCXS      := $(patsubst $(TEX_DIR)/%.tex,$(DOCX_DIR)/%.docx,$(SOURCES))

LATEXMK    ?= latexmk
LATEX_ENG  ?= lualatex
PANDOC     ?= pandoc

LATEXMK_FLAGS := -$(LATEX_ENG) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=$(PDF_DIR)

.PHONY: all pdf docx clean distclean

all: pdf docx

pdf: $(PDFS)

docx: $(DOCXS)

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
