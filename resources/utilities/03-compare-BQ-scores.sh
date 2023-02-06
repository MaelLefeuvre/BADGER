#!/usr/bin/env bash

# usage: 03-compare-BQ-scores.sh <sample-rmdup.bam> <sample-rmdup-rescaled.bam>

BAM_COL=11
PYTHON_CMD=$(cat << EOM
import sys
from math import floor
PHRED_SHIFT=33
for line in sys.stdin:
    line    = line.strip('\n').split(" ")
    
    left_chr = line[0]
    left_pos = int(line[1])
    left_nts = line[2]
    left_qlt = line[3]

    right_chr = line[4]
    right_pos = int(line[5])
    right_nts = line[6]
    right_qlt = line[7]

    if (left_pos != right_pos) or (left_chr != right_chr):
        print("Error: read nucleotide positions differ", file=sys.stderr)
        raise RuntimeError
        sys.exit(1)

    if (left_nts != right_nts):
        print("Error: read nucleotide sequence differ", file=sys.stderr)
        raise RuntimeError
        sys.exit(1)


    cmp_len = min(len(left_qlt), len(right_qlt))

    for (i, (a, b)) in enumerate(zip(left_qlt, right_qlt)):
        if a != b:

            # Find the relative position from the closest extremity (5' or 3')
            relpos = floor(cmp_len/2) - floor(abs(i - (cmp_len/ 2)))

            # Pick the value of the rescaled nucleotide.
            left_nuc  = left_nts[i]
            right_nuc = right_nts[i]

            # Ensure both nucleotides are the same
            if left_nuc != right_nuc:
                print("Error: sampled nucleotide differs", file=sys.stderr)
                raise RuntimeError
                sys.exit(1)

            print(f"{left_chr}:{left_pos+i}", left_nuc, a, ord(a)-PHRED_SHIFT, b, ord(b)-PHRED_SHIFT, relpos)
EOM
)

paste <(samtools view --no-header $1 | awk -v col=$BAM_COL '{print $3, $4, $10, $col}') <(samtools view --no-header $2 | awk -vcol=$BAM_COL '{print $3, $4, $10, $col}') \
| awk '{OFS=" "}($4!=$8){print $1, $2, $3, $4, $5, $6, $7, $8}' \
| python -c "${PYTHON_CMD}" \
| awk '{before_scores+=$4; after_scores+=$6; avg_relpos+=$7; print}END{print "Average BQ-score before rescaling:", before_scores/NR; print "Average BQ-score after rescaling :", after_scores/NR; print "Average relative position:", avg_relpos/NR}'
