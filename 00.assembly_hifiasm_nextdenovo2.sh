# genome assembly using hifiasm
# HiFi-only assembly
mkdir -p duck_hifi; hifiasm -o duck_hifi/duck_hifi -t 48 duck_female_HiFi.min1000.fastq.gz 2> duck_hifi/duck_hifi.log; for gfa in duck_hifi/duck_hifi*.gfa; do awk '/^S/{print ">"$2; print $3}' "$gfa" | seqkit seq -w 100 > "${gfa%.gfa}.fasta"; done

# ONT-only assembly
mkdir -p duck_ont; hifiasm -o duck_ont/duck_ont -t 48 --ont duck_ONT.min20k.q10.fastq.gz 2> duck_ont/duck_ont.log; for gfa in duck_ont/duck_ont*.gfa; do awk '/^S/{print ">"$2; print $3}' "$gfa" | seqkit seq -w 100 > "${gfa%.gfa}.fasta"; done

# HiFi + ONT hybrid assembly
mkdir -p duck_hifi_ont_hybrid; hifiasm -o duck_hifi_ont_hybrid/duck_hifi_ont_hybrid -t 48 --ul duck_ONT.min20k.q10.fastq.gz duck_female_HiFi.min1000.fastq.gz 2> duck_hifi_ont_hybrid/duck_hifi_ont_hybrid.log; for gfa in duck_hifi_ont_hybrid/duck_hifi_ont_hybrid*.gfa; do awk '/^S/{print ">"$2; print $3}' "$gfa" | seqkit seq -w 100 > "${gfa%.gfa}.fasta"; done



# genome assembly using nextdenovo2
echo "
[General]
job_type = local
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 10
input_type = raw
read_type = ont
input_fofn = input.fofn
workdir = assembly_out
[correct_option]
read_cutoff = 5k
genome_size = 1.3g
pa_correction = 10
sort_options = -m 4g -t 10
minimap2_options_raw = -t 10
correction_options = -p 10
[assemble_option]
minimap2_options_cns = -t 10
nextgraph_options = -a 1
" > run.cfg


# Write ONT read file paths to input.fofn
echo "$ont" > input.fofn


# Run NextDenovo
nextDenovo run.cfg
