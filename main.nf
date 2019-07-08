#!/usr/bin/env nextflow

// Initialize empty variables
params.sra = false
params.bioproject = false
params.output_folder = false

def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/nf-sra <ARGUMENTS>
    
    Arguments:
      --sra                         Accession for SRA sample to download
      --bioproject                  Accession for BioProject group to download
      --output_folder               Folder to place downloaded files

    Options:
      --help                        Print this help message

    All downloaded files will be named <ACCESSION>.fastq.gz. Paired-end files will be interleaved.

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

assert params.sra || params.bioproject: "Please specify --sra or --bioproject"
assert params.output_folder: "Please specify --output_folder"

// Either download a single SRA, or an entire BioProject
if ( params.bioproject ){
  assert !params.sra: "--sra cannot be specified together with --bioproject"

    // Get all of the SRA accessions for a BioProject
    process getBioProjectAccessions {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    cpus 4
    memory "8 GB"
    errorStrategy "retry"
    publishDir "${params.output_folder}"

    input:
    val bioproject from params.bioproject
    
    output:
    file "${bioproject}.csv" 
    stdout accession_ch

    afterScript "rm *"

    """
#!/usr/bin/env python3

import json
import requests
import pandas as pd

def get_sub_key(obj, path):
    for k in path:
        assert isinstance(obj, dict), "%s is not a dict" % (obj)
        assert k in obj, "%s not found in %s" % (k, json.dumps(obj, indent=4))
        obj = obj[k]
    return obj

def check_expected_response(d):
    assert isinstance(d, dict) and "data" in d, "Did not find a dict with 'data': %s" % (json.dumps(d, indent=4))
    assert len(d["data"]) > 0, "Found zero items in list: %s" % (json.dumps(d, indent=4))

# Use the MGnify API to search for the BioProject
bioproject_list = requests.get("https://www.ebi.ac.uk/metagenomics/api/v1/studies?accession=%s" % ("${bioproject}")).json()

check_expected_response(bioproject_list)
assert len(bioproject_list["data"]) == 1, "Found more than one BioProject with that name"
bioproject = bioproject_list["data"][0]

# Get the link to the samples for this BioProject
sample_list = requests.get(
    get_sub_key(
        bioproject,
        ["relationships", "samples", "links", "related"]
    )
).json()

check_expected_response(sample_list)

def get_sra(biosample):
    run = requests.get(
        get_sub_key(
            biosample,
            ["relationships", "runs", "links", "related"]
        )
    ).json()
    check_expected_response(run)
    for i in run["data"]:
        yield i["id"]


sra_list = []

for biosample in sample_list["data"]:
    biosample_desc = {
        k: v
        for k, v in get_sub_key(biosample, ["attributes"]).items()
        if v is not None and (isinstance(v, str) or isinstance(v, float) or isinstance(v, int))
    }

    for sra_accession in get_sra(biosample):
        sra_list.append({
            "run": sra_accession,
            **biosample_desc
        })

sra_list = pd.DataFrame(sra_list)
sra_list.to_csv("${bioproject}.csv", index=None)

print("\\n".join(sra_list["run"].tolist()))

    """
    }

}
else {
    if ( params.sra ){
        assert !params.bioproject: "--sra cannot be specified together with --bioproject"

        Channel.from(params.sra).set{ accession_ch }

    }
}

// Download the FASTQ files
process downloadSraFastq {
    container "quay.io/fhcrc-microbiome/get_sra:v0.4"
    cpus 4
    memory "8 GB"
    // errorStrategy "retry"
    publishDir "${params.output_folder}"

    input:
    val accession from accession_ch.splitText()

    output:
    file "*.fastq.gz"

    afterScript "rm -rf *"

"""
# Set the accession string
accession=\$(echo \"\"\"${accession}\"\"\" | tr -d '\\n')

# Get each read
echo "Get the FASTQ files"
fastq-dump --split-files --defline-seq '@\$ac.\$si.\$sg/\$ri' --defline-qual + --outdir \$PWD \$accession

r1=_1.fastq
r2=_2.fastq

# If there is a second read, interleave them
if [[ -s \$accession\$r2 ]]; then
    echo "Making paired reads"
    fastq_pair \$accession\$r1 \$accession\$r2
    
    echo "Interleave"
    paste <(gunzip -c \$accession\$r1.paired.fq) <(gunzip -c \$accession\$r2.paired.fq) | paste - - - - | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' | gzip -c > "\$accession.fastq.gz"
else
    echo "Compressing"
    mv \$accession\$r1 \$accession.fastq
    gzip \$accession.fastq
fi

rm -f \$accession\$r1 \$accession\$r2 \$accession\$r1.paired.fq \$accession\$r2.paired.fq

"""
}