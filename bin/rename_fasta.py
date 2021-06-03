#!/usr/bin/env python

import argparse
import csv
import sys
import re
import os


def _parse_name(seq):
    """Parse a fasta header and remove > and new lines.
    If [metadata is True] then parse the prophage metadata from
    the header.
    Current metadata: phage-circular and prophage-<start>:<end>
    """
    if not seq:
        return seq
    clean = seq.replace(">", "").replace("\n", "")

    clean = clean.replace("phage-circular", "")
    match = re.search(r"prophage-\d+:\d+", clean)
    prophage = match[0] if match else ""

    return clean.replace(prophage, "").strip(), \
        "phage-circular" if "phage-circular" in seq else "", prophage


def rename(args):
    """Rename a multi-fasta fasta entries with <name>.<counter> and store the
    mapping between new and old files in tsv (args.map)
    """
    print("Renaming " + args.input)
    with open(args.input, "r") as fasta_in:
        with open(args.output, "w") as fasta_out, open(args.map, "w") as map_tsv:
            count = 1
            tsv_map = csv.writer(map_tsv, delimiter="\t")
            tsv_map.writerow(["original", "renamed"])
            for line in fasta_in:
                if line.startswith(">"):
                    fasta_out.write(f">{args.prefix}{count}\n")
                    name, *_ = _parse_name(line)
                    tsv_map.writerow([name, f"{args.prefix}{count}"])
                    count += 1
                else:
                    fasta_out.write(line)
    print(f"Wrote {count} sequences to {args.output}.")

def main():
    """Multi fasta rename."""
    parser = argparse.ArgumentParser(
        description="Rename multi fastas.")
    parser.add_argument(
        "-i", "--input", help="indicate input FASTA file", required=False)
    parser.add_argument(
        "-m", "--map", help="Map current names with the renames", type=str,
        default="fasta_map.tsv")
    parser.add_argument(
        "-o", "--output", help="indicate output FASTA file", required=True)
    subparser = parser.add_subparsers()

    rename_parser = subparser.add_parser("rename")
    rename_parser.add_argument(
        "--prefix", help="string pre fasta count, i.e. default is seq such as seq1, seq2...",
        type=str, default="seq")
    rename_parser.set_defaults(func=rename)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
