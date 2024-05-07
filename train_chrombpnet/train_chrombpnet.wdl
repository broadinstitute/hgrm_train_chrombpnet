version development

workflow train_chrombpnet {
  input {
    File fragments
    File peaks
    File genome = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File chrom_sizes = "gs://joneslab-240402-village-training-data/chrombpnet_references/hg38.chrom.subset.sizes"
    File black_regions = "gs://joneslab-240402-village-training-data/chrombpnet_references/blacklist.bed.gz"
  }

  # TODO filter peaks against blacklist

  call get_background_regions {
    input:
    genome = genome,
    peaks = peaks,
    chrom_sizes = chrom_sizes,
    blacklist_regions = black_regions
  }

  call train_bias_model {
    input:
    fragments = fragments,
    genome = genome,
    peaks = peaks,
    chrom_sizes = chrom_sizes,
    peaks = peaks,
    non_peaks = get_background_regions.negatives,
    chr_folds = get_background_regions.folds
    }
}


##################################################
# GET BACKGROUND REGIONS
##################################################

task get_background_regions {
  input {
    File    genome
    File    peaks
    File    chrom_sizes
    File    blacklist_regions
  }

  Int disk_size = 20 + ceil(size(peaks, "GB"))

  command {
  mkdir -p output
  export TQDM_DISABLE=1  # clean up output
  chrombpnet prep splits -c ${chrom_sizes} -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op output/fold_0
  chrombpnet prep nonpeaks -g ${genome} -p ${peaks} -c ${chrom_sizes} -fl output/fold_0.json -br ${blacklist_regions} -o output/training
  ls -R output
  }

  output {
    File folds = "output/fold_0.json"
    File negatives = "output/training_negatives.bed"
  }

  runtime {
    docker: 'kundajelab/chrombpnet:latest'
    memory: "16 GB"
    bootDiskSizeGb: 20
    disks: "local-disk " + disk_size + " HDD"
  }
}


##################################################
# TRAIN BIAS MODEL
##################################################
task train_bias_model {
  input {
    File   fragments
    File   genome
    File   chrom_sizes
    File   peaks
    File   non_peaks
    File   chr_folds
  }

  Int disk_size = 100 + ceil(size(peaks, "GB")) + ceil(size(non_peaks, "GB"))

  command {
  set -euo pipefail
  mkdir -p output
  chrombpnet bias pipeline -ifrag ~{fragments} -d ATAC -g ~{genome} -c ~{chrom_sizes} -p ~{peaks} -n ~{non_peaks} -fl ~{chr_folds} -b 0.5 -o output/
  }

  output {
    File bias_model = "output/models/bias.h5"
    File log = "output/logs/bias.log"
    File evaluation = "output/evaluation/overall_report.pdf"
  }

  runtime {
    docker: 'kundajelab/chrombpnet:latest'
    memory: "256 GB"
    cpu: 24
    bootDiskSizeGb: 50
    disks: "local-disk " + disk_size + " HDD"
    gpuType: "nvidia-tesla-t4"
    gpuCount: 1
    nvidiaDriverVersion: "450.51.05"
    maxRetries: 1
  }
}
