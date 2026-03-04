#!/bin/bash
# Run the full entity resolution pipeline + script 17
# Usage: bash R/run_all_entity_resolution.sh

set -e
cd "/c/Users/Jared_Edgerton/Dropbox/Ecosystem_of_Aid"

RSCRIPT="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
PYTHON="/c/Users/Jared_Edgerton/AppData/Local/Programs/Python/Python312/python.exe"

echo "============================================"
echo "Step 1: Extract org names from networks"
echo "============================================"
"$RSCRIPT" R/entity_resolution_step1_extract_names.R

echo ""
echo "============================================"
echo "Step 2: Deterministic cleaning"
echo "============================================"
"$PYTHON" R/entity_resolution_step2_deterministic_clean.py

echo ""
echo "============================================"
echo "Step 3: Fuzzy deduplication"
echo "============================================"
"$PYTHON" R/entity_resolution_step3_fuzzy_dedup.py

echo ""
echo "============================================"
echo "Step 4a: Auto-review fuzzy candidates"
echo "============================================"
"$PYTHON" R/entity_resolution_step4a_review_fuzzy.py

echo ""
echo "============================================"
echo "Step 4c: Apply coder review (4 approved pairs only)"
echo "============================================"
"$PYTHON" R/entity_resolution_step4c_apply_coder_review.py

echo ""
echo "============================================"
echo "Step 4: Build crosswalk"
echo "============================================"
"$PYTHON" R/entity_resolution_step4_build_crosswalk.py

echo ""
echo "============================================"
echo "Step 5: Apply crosswalk to networks"
echo "============================================"
"$RSCRIPT" R/entity_resolution_step5_apply_crosswalk.R

echo ""
echo "============================================"
echo "Script 17: Build network + node panels"
echo "============================================"
"$RSCRIPT" R/17_network_topology.R

echo ""
echo "============================================"
echo "ALL DONE"
echo "============================================"
