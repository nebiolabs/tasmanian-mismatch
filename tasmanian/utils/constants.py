import numpy as np
class constants:
    SKIP_READS = {
        'indel':0,
        'mapquality':0,
        'fragLength':0,
        'softclips':0,
        'readlength':0,
        'others':0  # this could include cigar=* or mapq=255
    }
    SKIP_BASES = {
        'quality':0,
        'softclips':0
    }

    # set default parameters
    _UNMASK_GENOME=False  #don't unmask the genome, aka don't use lower letter from genome but rather keep them out in the error.log file'
    READ_LENGTH=1000
    PHRED = 20 + 33
    SOFTCLIP_BYPASS=0
    SKIP_INDEL=False
    MIN_LENGTH=0
    MAX_LENGTH=350
    TLEN=np.array([0,10000])
    ONT = False
    PWM  = False

    HELP = f"""
    required:
    --------
    -r|--reference-fasta

    optional:
    --------
    -u|--unmask-genome (convert masked bases to upper case and include them in the calculations - default=False)
    -q|--base-quality (default=20)
    -f|--filter-indel (exclude reads with indels default=False)
    -l|--filter-length (include only reads with x,y range of lengths, default={MIN_LENGTH},{MAX_LENGTH})
    -s|--soft-clip-bypass (Decide when softclipped base is correct(0). Don't use these bases(1). Force use them(2).  default=0)
    -m|--mapping-quality (minimum allowed mapping quality -defailt=0)
    -h|--help
    -g|--fragment-length (use fragments with these lengths ONLY)
    -o|--output-prefix (use this prefix for the output and logging files)
    -c|--confidence (number of bases in the confident region of the read) 
    -d|--debug (create a log file)
    -O|--ont (this is ONT data)
    -p|--picard-logic (normalize tables based on picard CollectSequencingArtifactMetrics logic)
    -P|--include-pwm
    """
