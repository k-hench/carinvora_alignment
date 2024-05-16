#!/usr/bin/awk -f

# Skip the header line
NR == 1 { print; next }

{
    # If this is the first data line or if the current line doesn't match the previous line,
    # print the previous line (if any) and store the current line as the new 'previous' line
    if (NR == 2 || $2 != prev_end || $4 != prev_val) {
        if (NR > 2) {
            print prev_chrom, prev_start, prev_end, prev_val
        }
        prev_chrom = $1
        prev_start = $2
        prev_end = $3
        prev_val = $4
    } else {
        # If the current line matches the previous line, just update the end position
        prev_end = $3
    }
}

# After the loop, print the last stored line
END {
    print prev_chrom, prev_start, prev_end, prev_val
}
