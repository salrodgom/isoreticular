# Reproducible-figure pipeline for the RHO isoreticular manuscript.
# `make all`            : rebuild every figure
# `make <target>`       : rebuild a single figure
# `make clean`          : wipe build/
# `make sync`           : copy build/*.pdf into the Overleaf images/ folder

GP        := gnuplot
PY        := python3

SCRIPTS_DIR := scripts
DATA_DIR    := data
BUILD_DIR   := build

OVERLEAF_IMG := /home/salvador/Dropbox/Aplicaciones/Overleaf/Structural\ changes\ in\ the\ RHO\ family\ of\ embedded\ isoreticular\ zeolites/images

# ---------- list of figures -----------------------------------------------
# Each <name> below is expected to map to scripts/<name>.gp producing
# build/<name>.eps, then converted to build/<name>.pdf.
FIGURES := \
    fig2_volume        fig2_cell \
    fig3a_volume       fig3b_delta       fig3e_softmode    fig3_phonons \
    fig5_elastic_G1    fig5_elastic_G2  fig5_elastic_G3 \
    fig5_kappa         fig5_isopressure \
    figS1_data_collapse \
    figS2_landau_wells_concept_a  figS2_landau_wells_concept_b \
    figS2_landau_wells_real_G1    figS2_landau_wells_real_G5 \
    figS3_landau_ueff \
    figS_landau_collapse_phi_eq  figS_landau_collapse_wells \
    figS4_g1_finite_size \
    figS5_g1_volume   figS5_g1_delta \
    figS5_g2_volume   figS5_g2_delta \
    figS6_density_energy \
    figS8_d8r_only_hist_G1  figS8_d8r_only_hist_G2  figS8_d8r_only_hist_G3 \
    figS8_d8r_only_hist_G4  figS8_d8r_only_hist_G5 \
    figS9_isolated_hist_G2  figS9_isolated_hist_G3 \
    figS9_isolated_hist_G4  figS9_isolated_hist_G5 \
    figS13_4mr_hist_G1  figS13_4mr_hist_G2  figS13_4mr_hist_G3 \
    figS13_4mr_hist_G4  figS13_4mr_hist_G5 \
    figS13_6mr_hist_G1  figS13_6mr_hist_G2  figS13_6mr_hist_G3 \
    figS13_6mr_hist_G4  figS13_6mr_hist_G5

EPSS := $(addprefix $(BUILD_DIR)/,$(addsuffix .eps,$(FIGURES)))

.PHONY: all clean sync help
all: $(EPSS)

# Generic rule: build/<name>.eps  <-  scripts/<name>.gp (+ any data/*.dat).
# The Overleaf pdflatex build ingests EPS directly via the LaTeX-side
# epstopdf package; no separate epstopdf binary is run here.
$(BUILD_DIR)/%.eps: $(SCRIPTS_DIR)/%.gp
	@mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && $(GP) ../$<

# Figure 3 data files come from the prep helper that reads the raw
# pressure-scan tables in initial_structures_RHO_isoreticular/.
$(BUILD_DIR)/fig3a_volume.eps $(BUILD_DIR)/fig3b_delta.eps: \
    $(DATA_DIR)/fig3_G1_volume_stable.dat

$(DATA_DIR)/fig3_G1_volume_stable.dat: $(SCRIPTS_DIR)/prep_fig3_volume_delta.py
	$(PY) $<

# Short aliases (so `make fig5_isopressure` works)
$(FIGURES): % : $(BUILD_DIR)/%.eps

clean:
	rm -rf $(BUILD_DIR)/*.eps $(BUILD_DIR)/*.pdf

sync: all
	./sync_to_overleaf.sh

help:
	@echo "Available figure targets:"
	@for f in $(FIGURES); do echo "  $$f"; done
