#!/usr/bin/env bash
# Single-pass splitter: writes tmp data while counting, then writes headers and gzips.
set -euo pipefail

COUNTS=GSE158055_covid19_counts.mtx.gz
BARC=GSE158055_covid19_barcodes.tsv.gz
FEAT=GSE158055_covid19_features.tsv.gz

for f in "$COUNTS" "$BARC" "$FEAT"; do [[ -f "$f" ]] || { echo "Missing $f"; exit 1; }; done

mkdir -p part1 part2

# sizes: <rows> <cols> <nnz>
read -r N_ROWS N_COLS NNZ < <(zcat "$COUNTS" | awk 'NR==1{next} /^%/{next} {print $1,$2,$3; exit}')
N1=${1:-$((N_COLS/2))}; (( N1>0 && N1<N_COLS )) || { echo "Bad N1=$N1"; exit 1; }
N2=$((N_COLS-N1))
echo "genes=$N_ROWS cells=$N_COLS nnz=$NNZ | part1=$N1 part2=$N2"

# split barcodes (no SIGPIPE) + copy features
zcat "$BARC" | awk -v n="$N1" 'NR<=n'    | gzip > part1/barcodes.tsv.gz
zcat "$BARC" | awk -v n="$N1" 'NR>n'     | gzip > part2/barcodes.tsv.gz
cp "$FEAT" part1/features.tsv.gz
cp "$FEAT" part2/features.tsv.gz

TMP1=$(mktemp part1/matrix.XXXXXX.tmp)
TMP2=$(mktemp part2/matrix.XXXXXX.tmp)

# SINGLE PASS: stream entries, write tmp data, count nnz, print progress
read -r NNZ1 NNZ2 < <(
  zcat "$COUNTS" | awk -v N1="$N1" -v T1="$TMP1" -v T2="$TMP2" '
    NR==1{next} /^%/{next} s==0{s=1; next}
    {
      r=$1; c=$2; v=$3; total++
      if (c<=N1) { print r, c, v >> T1; nnz1++ }
      else       { printf "%s %d %s\n", r, c-N1, v >> T2; nnz2++ }
      if (total % 10000000 == 0) {  # every 10M lines
        printf("[info] processed %d of ~%d (p1=%d p2=%d)\n", total, '"$NNZ"', nnz1, nnz2) > "/dev/stderr"
      }
    }
    END{ print (nnz1+0), (nnz2+0) }
  '
)

echo "nnz_part1=$NNZ1 nnz_part2=$NNZ2 (sum $((NNZ1+NNZ2)) / orig $NNZ)"

# write headers, then append gzipped data
{ echo "%%MatrixMarket matrix coordinate integer general"; echo "%"; echo "$N_ROWS $N1 $NNZ1"; } | gzip > part1/matrix.mtx.gz
{ echo "%%MatrixMarket matrix coordinate integer general"; echo "%"; echo "$N_ROWS $N2 $NNZ2"; } | gzip > part2/matrix.mtx.gz

gzip -c "$TMP1" >> part1/matrix.mtx.gz; rm -f "$TMP1"
gzip -c "$TMP2" >> part2/matrix.mtx.gz; rm -f "$TMP2"

echo "Done."
ls -lh part1 part2
