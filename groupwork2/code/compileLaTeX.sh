#!/bin/bash
# Compile LaTeX document

# Default to firstExample.tex if no argument is provided
if [ $# -eq 0 ]; then
    filename="firstExample"
else
    if [[ $1 == *.tex ]]; then
        filename=${1%.tex} # Remove file extension
    else
        filename=$1 # Use the provided name
    fi
fi

if [[ ! -f "$filename.tex" ]]; then
    echo "Error: File $filename.tex not found"
    exit 1
fi

pdflatex $filename.tex
bibtex $filename
pdflatex $filename.tex
pdflatex $filename.tex
evince $filename.pdf &

## Cleanup
rm *.aux
rm *.log
rm *.bbl
rm *.blg