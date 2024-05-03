version development

workflow train_chrombpnet {
  input {
    File fragments
    File peaks
    File? genome = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File? chrom_sizes = "gs://joneslab-240402-village-training-data/chrombpnet_references/hg38.chrom.subset.sizes"
    File? black_regions = "gs://joneslab-240402-village-training-data/chrombpnet_references/blacklist.bed.gz"
  }

  # TODO filter peaks against blacklist
  
  call get_background_regions {
    input: 
      genome = genome,
      peaks = peaks,
      chrom_sizes = chrom_sizes,
      blackRegions = black_regions
    }
}


##################################################
task get_background_regions {
  input {
    File    genome
    File    peaks
    File    chrom_sizes
    File    blackRegions
  }

  Int disk_size = 20 + size(peaks, “GB”)
  
  command {
  mkdir -p output
  chrombpnet prep splits -c ${chrom_sizes} -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op output/fold_0
  chrombpnet prep nonpeaks -g ${genome} -p ${peaks} -c ${chrom_sizes} -fl output/fold_0.json -br ${blackRegions} -o output/training
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
