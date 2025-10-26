#!/usr/bin/env python3、

#usage:  python3 align_seqs_fasta.py --seq1-from ****.fasta --seq2-from ****.fasta

from pathlib import Path
import argparse
import random
import difflib
from typing import Dict, Tuple, List

DEFAULT_FASTAS = ["407228326.fasta", "407228412.fasta", "E.coli.fasta"]

#read in fasta files
def read_fasta_long_short(filepath: str) -> Tuple[str, str, str]:
    """
      header (line 1), seq1 (concatenated long sequence), seq2 (last line).
    """
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with p.open("r") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()] #removes blank lines

    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Not a valid FASTA with a '>' header: {filepath}")
    if len(lines) < 3:
        raise ValueError(
            f"Expected at least 3 non-empty lines (header, long sequence, short sequence): {filepath}"
        )

    header = lines[0]
    seq1 = "".join(lines[1:-1]).upper()  # concatenate the long sequence pieces, and make sure they are in uppercases
    seq2 = lines[-1].upper()             # last line is the short sequence
    return header, seq1, seq2

#scoring alignment
def calculate_score(s1: str, s2: str, start: int) -> Tuple[int, str]:

    score = 0
    matched_chars: List[str] = []

    for i, base in enumerate(s2):
        j = i + start
        if j >= len(s1):
            break  # no further overlap
        if s1[j] == base:
            matched_chars.append("*") #match
            score += 1
        else:
            matched_chars.append("-") #mismatch

    return score, ("." * start) + "".join(matched_chars)

#find which is the best alignment by looping through all possible start positions in seq1
def best_alignment(seq1: str, seq2: str) -> Tuple[int, int, str, str, str]:
    """
    Find the best placement of seq2 (short) along seq1 (long).
    If multiple offsets tie, the last one found is kept.

    Returns:
      best_score, best_start, matched_line, aligned_s2_line, s1
    """
    best_score = -1
    best_start = 0
    best_matched_line = ""

    for start in range(len(seq1)):
        score, matched_line = calculate_score(seq1, seq2, start)
        if score > best_score:
            best_score = score
            best_start = start
            best_matched_line = matched_line

    aligned_s2_line = ("." * best_start) + seq2
    return best_score, best_start, best_matched_line, aligned_s2_line, seq1

#discover fasta files in working directory
def discover_fastas(candidates: List[str]) -> List[str]:
    """
    Return the subset of candidate filenames that exist in the working directory.
    """
    return [f for f in candidates if Path(f).exists()]

# fuzzy matching of file names, find closely matched files if no exact
def fuzzy_resolve(name: str, available: List[str]) -> str:
    if Path(name).exists():
        return name
    if not available:                                                       # Ensure file names match
        return name
    match = difflib.get_close_matches(name, available, n=1, cutoff=0.6)
    return match[0] if match else name

#which file is seq1 and which is seq2
def pick_pair(
    seq1_from: str | None, seq2_from: str | None, pool: List[str], rng: random.Random
) -> Tuple[str, str]:
    """
    Decide which file supplies seq1 and which supplies seq2.
    - If both provided -> use them (after fuzzy resolution).
    - If one provided -> pick the other at random from the pool.
    - If neither provided -> pick two at random, preferring different files.
    """
    available = discover_fastas(pool)
    if not available:
        raise FileNotFoundError(
            "No FASTA files found. Put them in the working directory or provide full paths."
        )

    # Fuzzy resolve names to available files (helps with minor typos)
    if seq1_from:
        seq1_from = fuzzy_resolve(seq1_from, available)
    if seq2_from:
        seq2_from = fuzzy_resolve(seq2_from, available)

    if seq1_from and seq2_from:
        return seq1_from, seq2_from

    if seq1_from and not seq2_from:
        choices = [f for f in available if f != seq1_from] or available
        return seq1_from, rng.choice(choices)

    if seq2_from and not seq1_from:
        choices = [f for f in available if f != seq2_from] or available
        return rng.choice(choices), seq2_from

    # neither provided -> draw two (prefer different files)
    if len(available) >= 2:
        a, b = rng.sample(available, 2)
        return a, b
    else:
        # Only one file present; both seq1 and seq2 come from it.
        return available[0], available[0]

# Main function to parse arguments and run alignment  
def main():
    parser = argparse.ArgumentParser(
        description="Align seq2 from one FASTA against seq1 from another."
    )
    parser.add_argument(
        "--seq1-from",
        dest="seq1_from",
        help="File to supply the long sequence (seq1).",
    )
    parser.add_argument(
        "--seq2-from",
        dest="seq2_from",
        help="File to supply the short sequence (seq2, last line).",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility."
    )
    parser.add_argument(
        "--pool",
        nargs="*",
        default=DEFAULT_FASTAS,
        help="Candidate FASTA files to consider when random choice is needed.",
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)

    # Decide which files provide seq1 and seq2
    file_seq1, file_seq2 = pick_pair(args.seq1_from, args.seq2_from, args.pool, rng)

    # Read both files and extract sequences
    header1, seq1_long, _ = read_fasta_long_short(file_seq1)
    header2, _, seq2_short = read_fasta_long_short(file_seq2)

    # Perform alignment (seq2 is aligned along seq1)
    best_score, best_start, matched_line, aligned_s2_line, s1 = best_alignment(
        seq1_long, seq2_short
    )

    # Report output
    print("\n=== Alignment summary ===")
    print(f"seq1 (long) from: {file_seq1}")
    print(f"  header ignored: {header1[1:]}")
    print(f"  length: {len(seq1_long)}")
    print(f"seq2 (short) from: {file_seq2}")
    print(f"  header ignored: {header2[1:]}")
    print(f"  length: {len(seq2_short)}")

    print("\nAlignment (matches='*', leading '.' shows the offset):")
    print(matched_line)
    print(aligned_s2_line)
    print(s1)
    print(f"\nBest score: {best_score} | Start position: {best_start}")


if __name__ == "__main__":
    main()