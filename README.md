# AllelicATAC

Usages:

Step 1. Please edit all the places marked XXXXXX in 0.mapping.sh. Make sure the paths and names of the FASTQ files are correct.

Step 2. ->If you want to run samples in batch, you can try using bsub.sh.
        ->If you don’t need batch mode, just run:
            bash 0.mapping.sh sample_name mat_label pat_label sample_label(where sample_label is the part of the filename before _R1_001.fastq.gz).

