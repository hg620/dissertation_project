#STEP 1: EXTRACTING SEQUENCES OF SELECTED GENES 
msa_file="input fasta file with sequences "
genes_presence_file="*** /scratch/gene_presence_counts.txt"
genes_boundaries_file="***/scratch/genes_boundaries.txt"
output_msa="*** /scratch/filtered_msa_100_lines.fasta"
temp_combined="*** /scratch/combined_sequences_random_100_lines.txt"

awk '$2 == 499' $genes_presence_file | shuf | head -n 1000 | cut -f1 | sort | uniq > selected_genes_random_100_lines.txt

# Extracting boundaries for selected genes
grep -Ff selected_genes_random_100_lines.txt $genes_boundaries_file > selected_genes_boundaries_random_100_lines.txt

# Initialising an empty associative array to store combined sequences
declare -A combined_sequences_100_lines

# Extracting and combine sequences for each genome
while read -r gene_name start end; do
    awk -v start="$start" -v end="$end" -v gene="$gene_name" -v OFS="" '
        BEGIN {header = ""}
        {
            if (NR % 2 == 1) {
                header = substr($0, 2)
            } else {
                seq = substr($0, start, end - start + 1)
                print ">" header "_" gene
                print seq
            }
        }
    ' $msa_file >> $temp_combined
done < selected_genes_boundaries_random_100_lines.txt

# Combining sequences per header in a dictionary, then write to an output file
#does not really work 
awk '
    /^>/ {header = $0; next}
    {
        seqs[header] = seqs[header] $0
    }
    END {
        for (header in seqs) {
            print header
            print seqs[header]
        }
    }
' $temp_combined > $output_msa

#STEP 2: CLEANUP AND CORRECTION. COMBINING GENES INTO ONE SEQUENCE PER SAMPLE 
input_file="*output of previous chunk of code"
declare -A sequences

#combining genes together to form an input for a tree 
while read -r line; do
    if [[ $line == \>* ]]; then
        header=$line
        sample=$(echo $header | cut -d' ' -f1)
        gene=$(echo $header | cut -d'_' -f2)
        combined_header=">${sample}"
    else
        sequences[$combined_header]+=$line
    fi
done < "$input_file"

for sample in "${!sequences[@]}"; do
    echo "$sample" >> "$output_file"
    echo "${sequences[$sample]}" >> "$output_file"
done
