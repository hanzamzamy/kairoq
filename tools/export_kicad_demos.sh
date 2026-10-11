#!/usr/bin/env bash
# Export the KiCad demonstration boards used as test boards into testdata/kicad/.
#
# Usage: tools/export_kicad_demos.sh [KICAD_SOURCE_DEMOS_DIR]
#
# Without an argument, the demos of the pinned KiCad release are downloaded into a
# temporary directory. Each board is exported with kicad-cli to Gerber X2 (with net
# attributes), Excellon and Gerber drill files, the Gerber job file, IPC-D-356,
# IPC-2581 (revision C), ODB++ and a pick-and-place CSV.
set -euo pipefail

KICAD_TAG=10.0.6
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=$ROOT/testdata/kicad

# name|path of the board below demos/|licence file of the demo directory, if it has one|
# IPC-2581 as plain XML ("xml") or zip-compressed ("zip")
BOARDS=(
    "ecc83|ecc83/ecc83-pp.kicad_pcb||xml"
    "pic_programmer|pic_programmer/pic_programmer.kicad_pcb||zip"
    "cm5_minima|cm5_minima/CM5_MINIMA_3.kicad_pcb|cm5_minima/LICENSE.txt|zip"
)

if [[ $# -ge 1 ]]; then
    DEMOS=$1
else
    TMP=$(mktemp -d)
    trap 'rm -rf "$TMP"' EXIT
    curl -sSL -o "$TMP/demos.tar.gz" \
        "https://gitlab.com/kicad/code/kicad/-/archive/$KICAD_TAG/kicad-$KICAD_TAG.tar.gz?path=demos"
    tar xzf "$TMP/demos.tar.gz" -C "$TMP"
    DEMOS=$TMP/kicad-$KICAD_TAG-demos/demos
fi

CLI_VERSION=$(kicad-cli version)
if [[ $CLI_VERSION != "$KICAD_TAG" ]]; then
    echo "warning: kicad-cli is $CLI_VERSION, the demos are from $KICAD_TAG" >&2
fi

for entry in "${BOARDS[@]}"; do
    IFS='|' read -r name pcb licence ipc <<<"$entry"
    src=$DEMOS/$pcb
    dst=$OUT/$name
    rm -rf "$dst"
    mkdir -p "$dst/gerber" "$dst/gerber-drill"
    # Copper, mask and outline only: other layers carry no information the tool uses.
    # Zones are not refilled, so that the Gerber and IPC-2581 exports share one fill.
    layers=$(sed -n 's/^\s*([0-9]\+ "\([A-Za-z0-9]*\.Cu\)".*/\1/p' "$src" | paste -sd, -)
    kicad-cli pcb export gerbers -o "$dst/gerber/" -l "$layers,F.Mask,B.Mask,Edge.Cuts" \
        "$src" >/dev/null
    kicad-cli pcb export drill -o "$dst/gerber/" --format excellon --excellon-separate-th \
        "$src" >/dev/null
    kicad-cli pcb export drill -o "$dst/gerber-drill/" --format gerber --excellon-separate-th \
        "$src" >/dev/null
    kicad-cli pcb export ipcd356 -o "$dst/gerber/$name.d356" "$src" >/dev/null
    kicad-cli pcb export pos -o "$dst/gerber/$name-pos.csv" --format csv --units mm "$src" >/dev/null
    if [[ $ipc == zip ]]; then
        kicad-cli pcb export ipc2581 -o "$dst/$name.ipc2581.zip" --compress --version C --units mm "$src" >/dev/null
    else
        kicad-cli pcb export ipc2581 -o "$dst/$name.xml" --version C --units mm "$src" >/dev/null
    fi
    # ODB++ with six decimals: the default of two rounds thicknesses and coordinates.
    kicad-cli pcb export odb -o "$dst/$name.odb.zip" --units mm --precision 6 "$src" >/dev/null
    if [[ -n $licence ]]; then
        cp "$DEMOS/$licence" "$dst/LICENSE.txt"
    fi
    {
        printf 'source: KiCad %s, demos/%s (https://gitlab.com/kicad/code/kicad)\n' "$KICAD_TAG" "$pcb"
        printf 'exported with: kicad-cli %s (tools/export_kicad_demos.sh)\n' "$CLI_VERSION"
        printf 'licence: CC BY-SA 4.0, the licence of the KiCad demo files (LICENSE.README of the KiCad source)\n'
        if [[ -n $licence ]]; then
            printf 'also: LICENSE.txt, the licence file of the demo directory\n'
        fi
        printf 'changes: exported to fabrication formats; no design changes\n'
    } >"$dst/SOURCE.txt"
    echo "$name: $(du -sh "$dst" | cut -f1)"
done
