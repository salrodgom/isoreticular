#!/usr/bin/env bash
# Unified PLUMED reweighting + FES generation pipeline.
# Usage:
#   bash run_reweight.sh <G_k>
# where <G_k> is 1, 2, 3, or 4. Must be run from the Multibaric directory
# containing MDNPT.{0,1,2}.xtc and COLVARS.{0,1,2}.
#
# For each target pressure P listed in the per-system PRESSURES_BAR array,
# generates fes_<P>.dat in cwd. Existing FES files at the same pressure are
# overwritten.

set -e

if [ $# -lt 1 ]; then
    echo "Usage: bash run_reweight.sh <G_k>  (1..4)"
    exit 1
fi
GK="$1"

# === Per-system parameters ===========================================
# Format: NSI, ENE_REF, VOL_REF, P_SIM_BAR, CELL_MIN, CELL_MAX, DELTA_MAX,
#          and PRESSURES_BAR = "p1 p2 p3 ..." (Fig 6 pressures in bar)
case "$GK" in
    1)  NSI=384;  ENE_REF=-4754832; VOL_REF=21.8; P_SIM_BAR=10000
        CELL_MIN=13.5; CELL_MAX=15.0; DELTA_MAX=2.5
        PRESSURES_BAR="10000 12000 12250 13000" ;;
    2)  NSI=270;  ENE_REF=-2971300; VOL_REF=14.3; P_SIM_BAR=7500
        CELL_MIN=23.5; CELL_MAX=25.0; DELTA_MAX=2.5
        PRESSURES_BAR="5000 7500 8200 9000" ;;
    3)  NSI=672;  ENE_REF=-8319800; VOL_REF=39.4; P_SIM_BAR=3500
        CELL_MIN=33.5; CELL_MAX=35.0; DELTA_MAX=3.0
        # Fig 6 has 1.0 GPa but OPES range is 1000-8000 bar so drop it
        PRESSURES_BAR="1000 5000 7000" ;;
    4)  NSI=1440; ENE_REF=-17830235; VOL_REF=86.0; P_SIM_BAR=1000
        CELL_MIN=43.5; CELL_MAX=44.5; DELTA_MAX=2.0
        # Fig 6 wants 2000-8000 bar but OPES range is 500-2000 only
        PRESSURES_BAR="500 1000 1500 2000" ;;
    *)  echo "G_$GK not supported. Use 1..4."; exit 1 ;;
esac

TEMP=300.0
KT=$(echo "scale=6; ${TEMP} * 0.008314410016255453" | bc -lq)
TPL_DIR="$(cd "$(dirname "$0")" && pwd)"
TPL="$TPL_DIR/template_plumed_REWEIGHT.tpl"

echo "=== Reweighting G_${GK} (N_Si=${NSI}, sim P=${P_SIM_BAR} bar) ==="
echo "Target pressures: $PRESSURES_BAR bar"
echo "Reference E=${ENE_REF} kJ/mol, V=${VOL_REF} nm^3"
echo ""

for P_BAR in $PRESSURES_BAR; do
    echo "--- P = ${P_BAR} bar ($(echo "scale=3; $P_BAR/10000" | bc) GPa) ---"
    P_SIM_PLU="0.06022140857*${P_SIM_BAR}"
    P_TGT_PLU="0.06022140857*${P_BAR}"

    # Build plumed_REWEIGHT files for each walker and run plumed driver
    for IDX in 0 1 2; do
        sed -e "s|__IDX__|${IDX}|g" \
            -e "s|__ENE_REF__|${ENE_REF}|g" \
            -e "s|__VOL_REF__|${VOL_REF}|g" \
            -e "s|__P_SIM__|${P_SIM_PLU}|g" \
            -e "s|__P_TARGET__|${P_TGT_PLU}|g" \
            -e "s|__TEMP__|${TEMP}|g" \
            "$TPL" > plumed_tmp.${IDX}
        plumed driver --plumed plumed_tmp.${IDX} --mf_xtc MDNPT.${IDX}.xtc --kt ${KT}
    done
    rm -f plumed_tmp.*

    # Concatenate the three walker outputs
    cat COLVAR_REWEIGHT.0 > COLVAR_cvs.${P_BAR}
    for IDX in 1 2; do
        sed '/#/d' COLVAR_REWEIGHT.${IDX} >> COLVAR_cvs.${P_BAR}
    done
    rm -f COLVAR_REWEIGHT.*

    # 2D FES on (cell, delta)
    python3 FES_from_reweighting.py \
        -f COLVAR_cvs.${P_BAR} -s 0.1,0.1 --cv cell,delta --kt ${KT} \
        --bias 7 --min ${CELL_MIN},0 --max ${CELL_MAX},${DELTA_MAX} --bin 150,150
    mv fes-rew.dat fes_${P_BAR}.dat
    rm -f bck.* COLVAR_cvs.${P_BAR}
    echo "  -> fes_${P_BAR}.dat written"
done

echo ""
echo "=== Done. Now run:"
echo "  python3 ${TPL_DIR}/extract_barriers.py > ${TPL_DIR}/barriers_extracted_G${GK}.txt"
echo "to compute ΔG* at each pressure."
