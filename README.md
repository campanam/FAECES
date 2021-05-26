# FAECES*: Fast and Accurate Enrichment of Canid Excrement for Species (\*and other analyses)  

This repository contains scripts to analyze noninvasive DNA from Mojave desert canids using the FAECES* assay.  

## Citation  
Please cite:  
Parker, L.D., Campana, M.G., Quinta, J.D., Cypher, B., Rivera, I., Fleischer, R.C., Ralls, K.S., Wilbert, T.R., Boarman, R., Boarman, W.I. & Maldonado, J.E. In prep. An efficient noninvasive method for simultaneous species, individual, and sex identification of sympatric Mojave Desert canids via in-solution SNP capture.  

## Usage  
`canid_sex.rb` performs sex determination using a previously called SNP VCF of the *ZFX/ZFY* region. Options include:  

    -i, --input [FILE]               Input VCF file
    -s, --species [FILE]             Species assignment file
    -k, --males [FILE]               List of known males to empirically calculate Y-allele frequency
        --mean                       Infer Y frequency using mean frequency across samples
    -z, --zscore [VALUE]             Infer mean Y frequency and set z-score bound
        --min                        Infer Y frequency using minimum frequency across samples
    -m, --minalleles [VALUE]         Minimum number of allele reads to call sex (Default = 10)
    -y, --yfreq [VALUE]              Expected Y-allele frequency (Default = 0.5)
    -a, --alpha [VALUE]              Alpha value (Default = 0.05)
        --yx [VALUE]                 Use specified Y-X allele ratio to determine sex rather than statistical test.
    -Z [VALUE]                       Use z-score cut-off to determine Y-X ratio.
    -h, --help                       Show help
