#!/usr/bin/env nextflow

// Initialize empty variables
params.sra = false
params.bioproject = false
params.output_folder = false
params.email = "sminot@fredhutch.org"

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
    container "quay.io/fhcrc-microbiome/biopython-pandas:latest"
    cpus 4
    memory "8 GB"
    errorStrategy "retry"
    publishDir "${params.output_folder}"

    input:
    val bioproject from params.bioproject
    val email from params.email
    
    output:
    file "${bioproject}.csv" 
    stdout accession_ch

    afterScript "rm *"

    """
#!/usr/bin/env python3

from Bio import Entrez
import pandas as pd
from time import sleep
import xmltodict
import json

Entrez.email = "${email}"

def get_sub_key(i, key_list):
    for k in key_list:
        assert isinstance(i, dict), "%s not a dict" % (json.dumps(i, indent=4))
        assert k in i, "%s not in %s" % (k, json.dumps(i, indent=4))
        i = i[k]
    return i

def get_samples_from_bioproject(bioproject):

    handle = Entrez.esearch(db="bioproject", term=bioproject)
    search_results = Entrez.read(handle)
    samples = Entrez.elink(dbfrom="bioproject", id=search_results["IdList"][0], linkname="bioproject_biosample")
    for element in Entrez.parse(samples):
        for linkset in get_sub_key(element, ["LinkSetDb"]):
            for sample in get_sub_key(linkset, ["Link"]):
                yield sample["Id"]
            
def get_sample_info(biosample_id):
    try:
        biosample = Entrez.efetch(db="biosample", id=biosample_id)
    except:
        sleep(5)
        biosample = Entrez.efetch(db="biosample", id=biosample_id)
    biosample = xmltodict.parse("".join([line for line in biosample]))
    
    d = {}
    
    for i in get_sub_key(biosample, ["BioSampleSet", "BioSample", "Attributes", "Attribute"]):
        if "@attribute_name" in i and "#text" in i:
            d[i["@attribute_name"]] = i["#text"]
    
    for i in get_sub_key(biosample, ["BioSampleSet", "BioSample", "Ids", "Id"]):
        for k in ["@db", "@db_label"]:
            if k in i and "#text" in i:
                d[i[k]] = i["#text"]
                
    for sra_accession in get_biosample_runs(biosample_id):
        d["Run"] = sra_accession
        yield d
        
def get_biosample_runs(biosample_id):
    try:
        handle = Entrez.elink(dbfrom="biosample", id=biosample_id, linkname="biosample_sra")
    except:
        sleep(5)
        handle = Entrez.elink(dbfrom="biosample", id=biosample_id, linkname="biosample_sra")
    search_results = xmltodict.parse("".join([line for line in handle]))
    
    links = get_sub_key(search_results, ["eLinkResult", "LinkSet", "LinkSetDb"])
    if isinstance(links, dict):
        links = [links]
    
    for link in links:
        try:
            run = Entrez.efetch(db="sra", id=get_sub_key(link, ["Link", "Id"]))
        except:
            sleep(5)
            run = Entrez.efetch(db="sra", id=get_sub_key(link, ["Link", "Id"]))
        run = xmltodict.parse("".join([line for line in run]))
        yield get_sub_key(run, ["EXPERIMENT_PACKAGE_SET", "EXPERIMENT_PACKAGE", "RUN_SET", "RUN", "@accession"])
    
df = []
for sample_id in get_samples_from_bioproject("${bioproject}"):
    for run in get_sample_info(sample_id):
        df.append(run)

df = pd.DataFrame(df)

df.to_csv("${bioproject}.csv", index=None)

print("\\n".join(df["Run"].tolist()))

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
    errorStrategy "retry"
    publishDir "${params.output_folder}"

    input:
    val accession from accession_ch.splitText().map{acc -> acc.trim()}

    output:
    file "*.fastq.gz"

    afterScript "rm -rf *"

"""
# Cache to the local folder
mkdir -p ~/.ncbi
mkdir cache
echo '/repository/user/main/public/root = "\$PWD/cache"' > ~/.ncbi/user-settings.mkfg

# Get each read
echo "Get the FASTQ files"
fastq-dump --split-files --defline-seq '@\$ac.\$si.\$sg/\$ri' --defline-qual + --outdir \$PWD ${accession}

r1=_1.fastq
r2=_2.fastq

# If there is a second read, interleave them
if [[ -s ${accession}\$r2 ]]; then
    echo "Making paired reads"
    fastq_pair ${accession}\$r1 ${accession}\$r2
    
    echo "Interleave"
    paste <(gunzip -c ${accession}\$r1.paired.fq) <(gunzip -c ${accession}\$r2.paired.fq) | paste - - - - | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' | gzip -c > "${accession}.fastq.gz"
else
    echo "Compressing"
    mv ${accession}\$r1 ${accession}.fastq
    gzip ${accession}.fastq
fi

rm -f ${accession}\$r1 ${accession}\$r2 ${accession}\$r1.paired.fq ${accession}\$r2.paired.fq

"""
}