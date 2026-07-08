#!/usr/bin/env bash
# Compile main6.tex and main10.tex, then build a latexdiff between them.
set -u
cd "$(dirname "$0")"

PDFLATEX=(pdflatex -interaction=nonstopmode -halt-on-error)

compile() {
    local base="$1"
    echo "==> Compiling ${base}.tex"
    "${PDFLATEX[@]}" "${base}.tex"  > "${base}.compile.log" 2>&1
    bibtex "${base}"               >> "${base}.compile.log" 2>&1
    "${PDFLATEX[@]}" "${base}.tex" >> "${base}.compile.log" 2>&1
    "${PDFLATEX[@]}" "${base}.tex" >> "${base}.compile.log" 2>&1
    local status=$?

    if [ "$status" -eq 0 ] && [ -f "${base}.pdf" ]; then
        echo "    OK -> ${base}.pdf"
        rm -f "${base}.aux" "${base}.out" "${base}.toc" "${base}.bbl" "${base}.blg" "${base}.compile.log"
    else
        echo "    FAILED - see ${base}.compile.log"
    fi
}

compile main6
compile main10

echo "==> Running latexdiff main6.tex main10.tex"
latexdiff main6.tex main10.tex > diff_main6_main10.tex 2> diff_main6_main10.latexdiff.log
if [ -s diff_main6_main10.latexdiff.log ]; then
    echo "    (latexdiff warnings logged in diff_main6_main10.latexdiff.log)"
else
    rm -f diff_main6_main10.latexdiff.log
fi

# latexdiff sometimes wraps table rows that are already commented out in
# main6.tex with \DIFdelbeginFL/\DIFdelendFL. Those markers' invisible \let
# side effects break booktabs' \noalign handling when placed directly next
# to \toprule/\midrule/\hline/\bottomrule. Strip them - no visible content
# is lost since the wrapped material is already invisible comments.
sed -i \
    -e 's/\\DIFdelbeginFL %DIFDELCMD/%DIFDELCMD/' \
    -e 's/\\DIFdelendFL \\toprule/\\toprule/' \
    -e 's/\\DIFdelendFL \\midrule/\\midrule/' \
    -e 's/\\DIFdelendFL \\hline/\\hline/' \
    -e 's/\\DIFdelendFL \\bottomrule/\\bottomrule/' \
    diff_main6_main10.tex

compile diff_main6_main10

echo "==> Done."
