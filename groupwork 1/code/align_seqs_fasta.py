#!/usr/bin/env python3

#usage:  python3 align_seqs_fasta.py ****.fasta ****.fasta
import os
import sys
import csv
from pathlib import Path
from typing import Tuple, List

# get script directory
script_dir = Path(__file__).parent
data_dir = script_dir / ".." / "data"
results_dir = script_dir / ".." / "results"


def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    seq = ''
    for line in lines:
        if not line.startswith('>'):
            seq += line.strip()
    return seq.upper()


def write_results_csv(filepath, file_seq1, file_seq2, best_score, best_start, best_matched_line, aligned_s2_line):
    """write alignment results to CSV file"""
    # convert to relative paths if they are absolute
    try:
        rel_file1 = Path(file_seq1).relative_to(Path.cwd())
    except ValueError:
        rel_file1 = Path(file_seq1).name
    
    try:
        rel_file2 = Path(file_seq2).relative_to(Path.cwd())
    except ValueError:
        rel_file2 = Path(file_seq2).name
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["file 1 is:", str(rel_file1)])
        writer.writerow([])  # empty row
        writer.writerow(["file 2 is:", str(rel_file2)])
        writer.writerow([])  # empty row
        writer.writerow(["best score is:", best_score])
        writer.writerow([])  # empty row
        writer.writerow(["start position is:", best_start])
        writer.writerow([])  # empty row
        writer.writerow(["best matched line is:", best_matched_line])
        writer.writerow([])  # empty row
        writer.writerow(["aligned s2 line is:", aligned_s2_line])


def default_arguments():
    
    # define default files using relative paths from script location
    default_file1 = data_dir / "407228412.fasta"
    default_file2 = data_dir / "407228326.fasta"
    
    if len(sys.argv) == 3:
        file_seq1 = sys.argv[1]
        file_seq2 = sys.argv[2]
    elif len(sys.argv) == 1:
        # use default files
        file_seq1 = str(default_file1.resolve())
        file_seq2 = str(default_file2.resolve())
    else:
        print("Usage: python3 align_seqs_fasta.py [seq1.fasta seq2.fasta]")
        print(f"If no files provided, defaults to {default_file1} and {default_file2}")
        sys.exit(1)
    
    # check if files exist
    if not Path(file_seq1).exists():
        print(f"Error: File not found: {file_seq1}")
        sys.exit(1)
    if not Path(file_seq2).exists():
        print(f"Error: File not found: {file_seq2}")
        sys.exit(1)
    
    # read sequences using read_fasta
    s1 = read_fasta(file_seq1)
    s2 = read_fasta(file_seq2)
    
    return file_seq1, file_seq2, s1, s2


# scoring alignment
def calculate_score(s1: str, s2: str, start: int) -> Tuple[int, str]:
    
    # ensure s1 is longer than s2
    l1 = len(s1)
    l2 = len(s2)
    if l1 >= l2:
        s1 = s1
        s2 = s2
    else:
        s1 = s2
        s2 = s1
        l1, l2 = l2, l1  # swap the two lengths

    score = 0
    matched_chars: List[str] = []

    for i, base in enumerate(s2):
        j = i + start
        if j >= l1:
            break  # no further overlap
        
        # check match and update both score and matched_chars
        is_match = (s1[j] == base)
        matched_chars.append("*" if is_match else "-")
        score += is_match

    return score, ("." * start) + "".join(matched_chars)

#find which is the best alignment by looping through all possible start positions in seq1
def best_alignment(seq1: str, seq2: str) -> Tuple[int, int, str, str]:
    """
    find the best placement of seq2 (short) along seq1 (long).
    if multiple positions have the same best score, keep the first one.

    returns:
      best_score: maximum number of matches
      best_start: offset position with best score
      best_matched_line: match pattern with '*' (match) and '-' (mismatch)
      aligned_s2_line: seq2 with leading dots showing offset
    """
    # find best alignment using max() - keeps first occurrence of best score
    find_best= max(
        (calculate_score(seq1, seq2, start) + (start,) for start in range(len(seq1))),
        key=lambda x: x[0]  # compare by score (first element)
    )
    best_score, best_matched_line, best_start = find_best

    aligned_s2_line = ("." * best_start) + seq2
    print(best_score, "\n", best_start, "\n", best_matched_line, "\n", aligned_s2_line)
    return best_score, best_start, best_matched_line, aligned_s2_line



# run alignment and print results
def run_alignment(file_seq1: str, file_seq2: str, s1: str, s2: str):
    
    best_score, best_start, matched_line, aligned_s2_line = best_alignment(s1, s2)
    
    print("\n=== Alignment summary ===")

    print("\nAlignment (matches='*', leading '.' shows the offset):")

    print("Best matched line is:", matched_line)
    print("Aligned seq2 line is:", aligned_s2_line)
    print(f"\nBest score: {best_score} | Start position: {best_start}")

    # write CSV
    os.makedirs(results_dir, exist_ok=True)
    csv_path = results_dir / "Alignment_results.csv"
    write_results_csv(csv_path, file_seq1, file_seq2, best_score, best_start, matched_line, aligned_s2_line)
    print(f"\nResults saved to: {csv_path}")  

# main function
def main():
    file_seq1, file_seq2, s1, s2 = default_arguments()
    run_alignment(file_seq1, file_seq2, s1, s2)


if __name__ == "__main__":
    main()