#!/usr/bin/env python3

# usage: python3 align_seqs_fasta.py --seq1-from ****.fasta --seq2-from ****.fasta

from pathlib import Path
import argparse
import random
import difflib

DEFAULT_FASTAS = ["407228326.fasta", "407228412.fasta", "E.coli.fasta"]

#shortened here, by merching (if not lines or not lines[0] and dlen(lines) < 3)
# and also removed (with p.open)
def read_fasta(filepath):
    """Return header, long sequence (all except last), short sequence (last line)."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">") or len(lines) < 3:
        raise ValueError(f"Invalid FASTA format: {filepath}")

    header = lines[0]
    seq1 = "".join(lines[1:-1]).upper()
    seq2 = lines[-1].upper()
    return header, seq1, seq2

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

def best_alignment(seq1, seq2):
    """Return the best alignment of seq2 along seq1."""
    best_score = -1
    for start in range(len(seq1)):
        score, matched_line = calculate_score(seq1, seq2, start)
        if score >= best_score: #merched here as well
            best_score, best_start, best_matched_line = score, start, matched_line

    aligned_s2 = "." * best_start + seq2
    return best_score, best_start, best_matched_line, aligned_s2, seq1

def discover_fastas(pool):
    """Return existing files from the pool."""
    return [f for f in pool if Path(f).exists()]

def fuzzy_resolve(name, available):
    """Return exact or closest filename match."""
    if Path(name).exists() or not available:
        return name
    matches = difflib.get_close_matches(name, available, n=1, cutoff=0.6)
    return matches[0] if matches else name

def pick_pair(seq1_from, seq2_from, pool, rng):
    """Decide which files provide seq1 and seq2."""
    available = discover_fastas(pool)
    if not available:
        raise FileNotFoundError("No FASTA files found in the working directory.")
    
    if seq1_from:
        seq1_from = fuzzy_resolve(seq1_from, available)
    if seq2_from:
        seq2_from = fuzzy_resolve(seq2_from, available)

    if seq1_from and seq2_from:
        return seq1_from, seq2_from
    if seq1_from:
        choices = [f for f in available if f != seq1_from] or available
        return seq1_from, rng.choice(choices)
    if seq2_from:
        choices = [f for f in available if f != seq2_from] or available
        return rng.choice(choices), seq2_from
    return tuple(rng.sample(available, 2)) if len(available) >= 2 else (available[0], available[0])

def main():
    parser = argparse.ArgumentParser(description="Align seq2 from one FASTA against seq1 from another.")
    parser.add_argument("--seq1-from", help="File for the long sequence (seq1).")
    parser.add_argument("--seq2-from", help="File for the short sequence (seq2).")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility.")
    parser.add_argument("--pool", nargs="*", default=DEFAULT_FASTAS, help="Candidate FASTA files.")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    file_seq1, file_seq2 = pick_pair(args.seq1_from, args.seq2_from, args.pool, rng)

    header1, seq1_long, _ = read_fasta(file_seq1)
    header2, _, seq2_short = read_fasta(file_seq2)

    best_score, best_start, matched_line, aligned_s2_line, seq1 = best_alignment(seq1_long, seq2_short)

    print("\n=== Alignment summary ===")
    print(f"seq1 (long) from: {file_seq1} | length: {len(seq1_long)} | header ignored: {header1[1:]}")
    print(f"seq2 (short) from: {file_seq2} | length: {len(seq2_short)} | header ignored: {header2[1:]}")
    print("\nAlignment (matches='*', leading '.' shows the offset):")
    print(matched_line)
    print(aligned_s2_line)
    print(seq1)
    print(f"\nBest score: {best_score} | Start position: {best_start}")

if __name__ == "__main__":
    main()
