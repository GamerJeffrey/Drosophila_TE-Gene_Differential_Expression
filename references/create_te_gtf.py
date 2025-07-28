#!/usr/bin/env python3
import sys

def create_te_gtf(rm_file, gtf_file):
    with open(rm_file, 'r') as f, open(gtf_file, 'w') as out:
        for i, line in enumerate(f):
            if i < 3:
                continue
            fields = line.strip().split()
            if len(fields) < 15:
                continue
            chrom, start, end = fields[4], fields[5], fields[6]
            strand = '+' if fields[8] == '+' else '-'
            te_name, te_class = fields[9], fields[10]
            te_id = f"{te_name}_{chrom}_{start}_{end}_{strand}"
            attrs = f'gene_id "{te_id}"; transcript_id "{te_id}"; gene_name "{te_name}"; te_class "{te_class}";'
            out.write(f"{chrom}\tRepeatMasker\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <repeatmasker.out> <output.gtf>", file=sys.stderr)
        sys.exit(1)
    create_te_gtf(sys.argv[1], sys.argv[2])
