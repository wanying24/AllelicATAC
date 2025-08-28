#!/bin/bash

DATA='XXXXX' ####datadir for fastq

declare -A key_sample
declare -A key_pat
make_unique() {
    echo "$1" | tr ' ' '\n' | sort -u | tr '\n' ' '
}
for file in ${DATA}CC*_R1_001.fastq.gz; do
    filename=$(basename "$file")
    PAT=${filename%_R*}
    IFS='_' read -ra PARTS <<< "$filename"
    key="${PARTS[0]//J/}"
    key_pat[${PARTS[1]}]=$PAT
    key="${key//x/_}"
    if [ -n "${PARTS[1]}" ]; then
        key_sample[$key]+="${PARTS[1]} "
        
    fi
done

for s in "${!key_sample[@]}"; do
    SAMPLES=$(make_unique "${key_sample[$s]}")
    part1="${s%%_*}"
    part2="${s##*_}"
    for e in ${SAMPLES};do 
        CC=${key_pat[$e]}
        script_path="${e}_job.sh"
        cat > "$script_path" <<EOF
#!/bin/bash
#BSUB -q XXXXX
#BSUB -G XXXXX
#BSUB -g XXXXX
#BSUB -a 'docker(wuw2024/washu:atac0708)'
#BSUB -n 16
#BSUB -M 100GB
#BSUB -R 'select[mem>32GB && tmp>32GB] rusage[mem=60GB, tmp=60GB] span[hosts=1]'
#BSUB -W 36000
#BSUB -e /storage1/fs1/yeli/Active/wanying/log/%J_stderr.log
export PATH="/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/opt/conda/envs/nf-core-atacseq-1.2.2/bin/:\$PATH"
cd XXXXXXXXXXX
bash 0.mapping.sh ${e} ${part1} ${part2} ${CC}
EOF
        echo "Generated job script for sample $e at $script_path"
    done
done
