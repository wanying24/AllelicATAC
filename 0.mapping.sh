#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
trap 'echo "[ERROR] at line $LINENO"; exit 1' ERR

# Usage:
#   ./0.mapping.sh SAMPLE_NAME MAT_LABEL PAT_LABEL
# Example:
#   ./0.mapping.sh JJ320 CC001 CC074 JJ320

NAME=${1:?need NAME}
CC1=${2:?need CC1}
CC2=${3:?need CC2}
PAT=${4:?need PAT}


###Add by yourself
DATA='XXXX' ###fastq_dir
OUT='XXXX' ###out_dir
GENOME='XXXX' ###genome_dir
SCRIPTS='XXXXX' ###script_dir
T=24                             # threads
PICARD_JAVA_OPTS="-Xmx32g"      

F1="${DATA}/${PAT}_R1_001.fastq.gz" ###double check fastq name
F2="${DATA}/${PAT}_R2_001.fastq.gz"

# ---------- helpers ----------
mkdir_p() { [[ -d "$1" ]] || { mkdir -p "$1"; echo "[MKDIR] $1"; }; }
nreads() { zcat "$1" | awk 'NR%4==1{c++} END{print c+0}'; }  

# ---------- dirs ----------
mkdir_p "${OUT}/1.clean"
mkdir_p "${OUT}/2.statics"
mkdir_p "${OUT}/3.pesduoBam"

cd "${OUT}/3.pesduoBam"

STATS="${OUT}/2.statics/${NAME}_unique_file.txt"
echo -e "raw reads\t$(nreads "${F1}")" > "${STATS}"

# ---------- 1) Trim Galore ----------
# Nextera  'CTGTCTCTTATACACATCT'   #####double check
trim_galore -j "${T}" -q 36 --length 36 --max_n 3 --stringency 3 --phred33 \
  --paired -a 'CTGTCTCTTATACACATCT' --fastqc \
  -o "${OUT}/1.clean" "${F1}" "${F2}"

R1_CLEAN=$(ls "${OUT}/1.clean/" | grep -E "${PAT}.*_val_1\.f(ast)?q\.gz$" | head -n1)
R2_CLEAN=$(ls "${OUT}/1.clean/" | grep -E "${PAT}.*_val_2\.f(ast)?q\.gz$" | head -n1)
[[ -n "${R1_CLEAN:-}" && -n "${R2_CLEAN:-}" ]] || { echo "[ERR] trimmed files not found"; exit 1; }

mv -f "${OUT}/1.clean/${R1_CLEAN}" "${OUT}/1.clean/${NAME}_R1.fastq.gz"
mv -f "${OUT}/1.clean/${R2_CLEAN}" "${OUT}/1.clean/${NAME}_R2.fastq.gz"

F1C="${OUT}/1.clean/${NAME}_R1.fastq.gz"
F2C="${OUT}/1.clean/${NAME}_R2.fastq.gz"
echo -e "after QC (clean) reads\t$(nreads "${F1C}")" >> "${STATS}"

# ---------- 2) Bowtie2 index ----------
PSEUDO="${GENOME}/${CC1}x${CC2}_genome.fa"
if ls "${PSEUDO}".*.bt2* 1>/dev/null 2>&1; then
  echo "[INFO] Bowtie2 index exists."
else
  echo "[INFO] Building Bowtie2 index ..."
  bowtie2-build --threads "${T}" "${PSEUDO}" "${PSEUDO}"
fi

# ---------- 3) Align ----------
BAM_RAW="${NAME}_${CC1}x${CC2}_genome.bam"
BAM_SORT="${NAME}_${CC1}x${CC2}_genome_sorted.bam"
BAM_DEDUP="${NAME}_${CC1}x${CC2}_genome_sorted.dedup.bam"
BAM_NAME_SORT="${NAME}_${CC1}x${CC2}_genome_srt.bam"

echo "[INFO] Aligning ..."

bowtie2 --very-sensitive -X 2000 -k 10 -p "${T}" \
  -x "${PSEUDO}" -1 "${F1C}" -2 "${F2C}" \
| samtools view -@ "${T}" -bS - > "${BAM_RAW}"

echo "[INFO] Sorting ..."
samtools sort -@ "${T}" -o "${BAM_SORT}" "${BAM_RAW}"
samtools index -@ "${T}" "${BAM_SORT}"


samtools flagstat -@ "${T}" "${BAM_SORT}" > "${OUT}/2.statics/${NAME}_${CC1}x${CC2}_flagstat.txt"
samtools idxstats "${BAM_SORT}" > "${OUT}/2.statics/${NAME}_${CC1}x${CC2}_idxstats.txt"
samtools stats -@ "${T}" "${BAM_SORT}" > "${OUT}/2.statics/${NAME}_${CC1}x${CC2}_samstats.txt"

# ---------- 4) MarkDuplicates ----------
echo "[INFO] MarkDuplicates ..."
picard ${PICARD_JAVA_OPTS} MarkDuplicates \
  -I "${BAM_SORT}" \
  -O "${BAM_DEDUP}" \
  -M "${OUT}/2.statics/${NAME}_${CC1}x${CC2}_genome.dedup_metrics" \
  --REMOVE_DUPLICATES true

samtools index -@ "${T}" "${BAM_DEDUP}"
samtools flagstat -@ "${T}" "${BAM_DEDUP}" > "${OUT}/2.statics/${NAME}_${CC1}x${CC2}_flagstat.dedup.txt"

# ---------- 5) Name-sort for splitting ----------
echo "[INFO] Name-sorting ..."
samtools sort -@ "${T}" -n -o "${BAM_NAME_SORT}" "${BAM_DEDUP}"

# ---------- 6) Split by allele ----------
echo "[INFO] Splitting by allele ..."
python "${SCRIPTS}/1.split_reads.py" \
  "${BAM_NAME_SORT}" "${NAME}_${CC1}" "${NAME}_${CC2}"

echo "[DONE] ${NAME}  ${CC1}x${CC2}"
