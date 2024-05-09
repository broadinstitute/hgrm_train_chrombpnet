version development

workflow train_chrombpnet {
  input {
    String file_prefix
    File fragments
    File peaks
    File genome = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File chrom_sizes = "gs://joneslab-240402-village-training-data/chrombpnet_references/hg38.chrom.subset.sizes"
    File blacklist_regions = "gs://joneslab-240402-village-training-data/chrombpnet_references/blacklist.bed.gz"
    File? pretrained_bias_model
    Boolean run_training = true
  }

  call filter_peaks {
    input:
    peaks = peaks,
    chrom_sizes = chrom_sizes,
    blacklist_regions = blacklist_regions,
    file_prefix = file_prefix
  }

  call get_background_regions {
    input:
    genome = genome,
    peaks = filter_peaks.filtered_peaks,
    chrom_sizes = chrom_sizes,
    blacklist_regions = blacklist_regions,
    file_prefix = file_prefix
  }

  call generate_signal_bigwig {
    input:
    fragments = fragments,
    genome = genome,
    chrom_sizes = chrom_sizes,
    file_prefix = file_prefix
  }

  if (run_training) {
    if (! defined(pretrained_bias_model)) {
      call train_bias_model {
        input:
        fragments = fragments,
        genome = genome,
        peaks = filter_peaks.filtered_peaks,
        chrom_sizes = chrom_sizes,
        non_peaks = get_background_regions.negatives,
        chr_folds = get_background_regions.folds,
        file_prefix = file_prefix
      }
    }

    call train_factorized_model {
      input:
      bias_model = select_first([pretrained_bias_model, train_bias_model.bias_model]),  # ick
      fragments = fragments,
      genome = genome,
      chrom_sizes = chrom_sizes,
      peaks = filter_peaks.filtered_peaks,
      non_peaks = get_background_regions.negatives,
      chr_folds = get_background_regions.folds,
      file_prefix = file_prefix
    }
  }

  output {
    File cbpn_filtered_peaks = filter_peaks.filtered_peaks
    File cbpn_background_peaks = get_background_regions.negatives
    File cbpn_signal_bigwig = generate_signal_bigwig.bigwig
    File? cbpn_bias_model = train_bias_model.bias_model
    File? cbpn_bias_model_evaluation = train_bias_model.evaluation
    File? cbpn_chrombpnet_model = train_factorized_model.chrombpnet_model
    File? cbpn_chrombpnet_model_nobias = train_factorized_model.chrombpnet_model_nobias
    File? cbpn_bias_model_scaled = train_factorized_model.bias_model_scaled
    File? cbpn_model_evaluation = train_factorized_model.evaluation
  }
}


##################################################
# FILTER PEAKS
##################################################

task filter_peaks {
  input {
    File    peaks
    File    chrom_sizes
    File    blacklist_regions
    String  file_prefix
  }

  Int disk_size = 20 + 3 * ceil(size(peaks, "GB"))

  command {
  mkdir -p outputs
  bedtools slop -i ~{blacklist_regions} -g ~{chrom_sizes} -b 1057 > temp.bed
  # filter peaks that overlap blacklist, and also trim to 10 columns (narrowPeak format)
  bedtools intersect -v -a ~{peaks} -b temp.bed  | cut -f 1-10 > outputs/~{file_prefix}.peaks_no_blacklist.bed
  }

  output {
    File filtered_peaks = "outputs/~{file_prefix}.peaks_no_blacklist.bed"
  }

  runtime {
    docker: "quay.io/biocontainers/bedtools:2.24--1"
    memory: "10 GB"
    bootDiskSizeGb: 20
    disks: "local-disk " + disk_size + " HDD"
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
    String  file_prefix
  }

  Int disk_size = 20 + ceil(size(peaks, "GB"))

  command {
  pip install -U tqdm
  mkdir -p output
  export TQDM_DISABLE=1  # clean up output
  chrombpnet prep splits -c ${chrom_sizes} -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op output/fold_0
  chrombpnet prep nonpeaks -g ${genome} -p ${peaks} -c ${chrom_sizes} -fl output/fold_0.json -br ${blacklist_regions} -o output/training
  cp output/training/training_negatives.bed output/training/${file_prefix}.training_negatives.bed
  }

  output {
    File folds = "output/fold_0.json"
    File negatives = "output/~{file_prefix}.training_negatives.bed"
  }

  runtime {
    docker: 'kundajelab/chrombpnet:latest'
    memory: "16 GB"
    bootDiskSizeGb: 20
    disks: "local-disk " + disk_size + " HDD"
  }
}


##################################################
# GENERATE SIGNAL BIGWIG
##################################################

task generate_signal_bigwig {
  input {
    File    fragments
    File    genome
    File    chrom_sizes
    String  file_prefix
  }

  Int disk_size = 20 + 3 * ceil(size(fragments, "GB"))

  command {
  mkdir -p output
  python /scratch/chrombpnet/chrombpnet/helpers/preprocessing/reads_to_bigwig.py \
         -ifrag ${fragments} -d ATAC -c ${chrom_sizes} -g ${genome} -op output/${file_prefix}_ATAC
  }

  output {
    File bigwig = "output/~{file_prefix}_ATAC_unstranded.bw"
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
    String file_prefix
  }

  Int disk_size = 100 + ceil(size(peaks, "GB")) + ceil(size(non_peaks, "GB"))

  command {
  set -euo pipefail
  pip install -U tqdm
  export TQDM_MININTERVAL=30
  mkdir -p output
  chrombpnet bias pipeline -ifrag ~{fragments} -d ATAC -g ~{genome} -c ~{chrom_sizes} -p ~{peaks} -n ~{non_peaks} -fl ~{chr_folds} -b 0.5 -o output/
  mv output/models/bias.h5 output/models/~{file_prefix}.bias.h5
  mv output/models/overall_report.pdf output/models/~{file_prefix}.overall_report.pdf
  }

  output {
    File bias_model = "output/models/~{file_prefix}.bias.h5"
    File log = "output/logs/bias.log"
    File evaluation = "output/evaluation/~{file_prefix}.overall_report.pdf"
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
  }
}


##################################################
# TRAIN FACTORIZED MODEL
##################################################
task train_factorized_model {
  input {
    File   bias_model
    File   fragments
    File   genome
    File   chrom_sizes
    File   peaks
    File   non_peaks
    File   chr_folds
    String   file_prefix
  }

  Int disk_size = 100 + ceil(size(peaks, "GB")) + ceil(size(non_peaks, "GB")) + ceil(size(bias_model, "GB"))

  command {
  set -euo pipefail
  mkdir -p output
  pip install -U tqdm
  export TQDM_MININTERVAL=30
  chrombpnet pipeline \
        -ifrag ~{fragments} \
        -d ATAC \
        -g ~{genome} \
        -c ~{chrom_sizes} \
        -p ~{peaks} \
        -n ~{non_peaks} \
        -fl ~{chr_folds} \
        -b ~{bias_model} \
        -o output/

  mv output/models/chrombpnet.h5               output/models/~{file_prefix}.chrombpnet.h5
  mv output/models/chrombpnet_nobias.h5        output/models/~{file_prefix}.chrombpnet_nobias.h5
  mv output/models/bias_model_scaled.h5        output/models/~{file_prefix}.bias_model_scaled.h5
  mv output/logs/chrombpnet.log                output/logs/~{file_prefix}.chrombpnet.log
  mv output/evaluation/overall_report.pdf      output/evaluation/~{file_prefix}.overall_report.pdf
  mv output/evaluation/chrombpnet_metrics.json output/evaluation/~{file_prefix}.chrombpnet_metrics.json
  }

  output {
    File chrombpnet_model =         "output/models/~{file_prefix}.chrombpnet.h5"
    File chrombpnet_model_nobias =  "output/models/~{file_prefix}.chrombpnet_nobias.h5"
    File bias_model_scaled =        "output/models/~{file_prefix}.bias_model_scaled.h5"
    File log =                      "output/logs/~{file_prefix}.chrombpnet.log"
    File evaluation =               "output/evaluation/~{file_prefix}.overall_report.pdf"
    File metrics =                  "output/evaluation/~{file_prefix}.chrombpnet_metrics.json"
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
   }
}
