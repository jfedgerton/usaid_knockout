###############################################################################
# Step 4c: Apply coder review — only accept the 4 approved fuzzy matches
#
# The coder reviewed all fuzzy candidates and approved only these 4 pairs:
#   1. ministere des travaux publics transports et communications
#      -> ministre des travaux publics transports et communications
#   2. international legal assisitance consortium
#      -> international legal assistance consortium
#   3. fond d assistance economique et sociale
#      -> fonds d assistance economique et sociale
#   4. ministere des travaux publics et de l entretien routier
#      -> ministre des travaux publics et de i entretien routier
#
# All other fuzzy candidates are rejected.
###############################################################################

import csv

INPUT = "Data/entity_resolution/reviewed_fuzzy_matches.csv"
OUTPUT = INPUT  # overwrite

# The 4 approved pairs (as tuples of name_a, name_b)
APPROVED_PAIRS = {
    ("ministere des travaux publics transports et communications",
     "ministre des travaux publics transports et communications"),
    ("international legal assisitance consortium",
     "international legal assistance consortium"),
    ("fond d assistance economique et sociale",
     "fonds d assistance economique et sociale"),
    ("ministere des travaux publics et de l entretien routier",
     "ministre des travaux publics et de i entretien routier"),
}

rows = []
with open(INPUT, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

n_accepted = 0
n_rejected = 0

for row in rows:
    pair = (row["name_a"].strip(), row["name_b"].strip())
    pair_rev = (row["name_b"].strip(), row["name_a"].strip())

    if pair in APPROVED_PAIRS or pair_rev in APPROVED_PAIRS:
        row["accept"] = "yes"
        row["review_note"] = "coder_approved"
        # Canonical: use the correct spelling (longer/standard form)
        # For typos, prefer the correctly-spelled version
        if "assisitance" in row["name_a"]:
            row["canonical_name"] = row["name_b"]  # "assistance" (correct)
        elif "ministere" in row["name_a"] and "ministre" in row["name_b"]:
            row["canonical_name"] = row["name_a"]  # "ministere" (correct)
        elif "fond " in row["name_a"] and "fonds " in row["name_b"]:
            row["canonical_name"] = row["name_b"]  # "fonds" (correct)
        n_accepted += 1
    else:
        row["accept"] = "no"
        if row["review_note"] not in ("country_name_swap", "differentiating_word",
                                       "word_count_mismatch", "short_names_low_score"):
            row["review_note"] = "coder_rejected"
        n_rejected += 1

with open(OUTPUT, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"Coder review applied:")
print(f"  Accepted: {n_accepted}")
print(f"  Rejected: {n_rejected}")
print(f"  Total: {len(rows)}")

# Verify the accepted ones
print(f"\nAccepted pairs:")
for row in rows:
    if row["accept"] == "yes":
        print(f"  {row['name_a']}")
        print(f"    -> {row['name_b']}")
        print(f"    canonical: {row['canonical_name']}")
