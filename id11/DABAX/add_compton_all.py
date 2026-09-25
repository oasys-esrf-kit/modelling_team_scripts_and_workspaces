#
# Extend add_compton.py's Si-only Compton correction to every element in
# f1f2_Windt.dat, writing the corrected f2 for all elements to
# f1f2_WindtWithCompton.dat
#
import re
import numpy
import xraylib


INPUT_FILE = "f1f2_Windt.dat"
OUTPUT_FILE = "f1f2_WindtWithCompton.dat"


def correct_f2(symbol, energies_eV, f2):
    """Multiply f2 by (photoelectric + compton) / photoelectric, per element.

    Same approach as the original Si-only script: energies outside the
    xraylib tables raise an exception, in which case the previous valid
    ratio is reused (this mainly affects the low-energy points below the
    element's absorption tables).
    """
    Z = xraylib.SymbolToAtomicNumber(symbol)

    energies_keV = energies_eV * 1e-3

    ratios = []
    ratio = 1.0
    for E in energies_keV:
        try:
            compton = xraylib.CS_Compt(Z, E)
            photoelectric = xraylib.CS_Photo(Z, E)
            ratio = (photoelectric + compton) / photoelectric if photoelectric != 0 else ratio
            ratios.append(ratio)
        except Exception:
            ratios.append(ratio)

    ratios = numpy.array(ratios)
    return f2 * ratios


if __name__ == "__main__":
    raw = open(INPUT_FILE).read()

    first_section_start = raw.find("#S ")
    preamble = raw[:first_section_start]

    # split into per-element blocks, each starting with a "#S <index> <Symbol>" line
    section_starts = [m.start() for m in re.finditer(r"(?m)^#S\s+\d+\s+\S+", raw)]
    section_starts.append(len(raw))

    out_chunks = [preamble]

    for i in range(len(section_starts) - 1):
        block = raw[section_starts[i]:section_starts[i + 1]]
        lines = block.splitlines()

        header_line = lines[0]
        symbol = header_line.split()[2]

        header_lines = []
        data_lines = []
        for line in lines[1:]:
            if line.startswith("#") or not line.strip():
                header_lines.append(line)
            else:
                data_lines.append(line)

        data = numpy.array([[float(v) for v in line.split()] for line in data_lines])
        energies_eV, f1, f2 = data[:, 0], data[:, 1], data[:, 2]

        f2_corrected = correct_f2(symbol, energies_eV, f2)

        print("Processing %3s (%d energy points)" % (symbol, energies_eV.size))

        out_lines = [header_line] + header_lines
        for e, a, b in zip(energies_eV, f1, f2_corrected):
            out_lines.append("%s %s %s" % (e, a, b))

        out_chunks.append("\n".join(out_lines) + "\n")

    with open(OUTPUT_FILE, "w") as f:
        f.write("".join(out_chunks))

    print("Written:", OUTPUT_FILE)
