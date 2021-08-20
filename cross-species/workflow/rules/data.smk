import pandas as pd

targets = (
	pd.read_csv(config["targets"], sep="\t", dtype={"species": str})
	.set_index("species", drop=False)
	.sort_index()
)

tsa_targets = targets.loc[lambda targets: targets['source'] == "tsa"]
ensembl_targets = targets.loc[lambda targets: targets['source'] == "ensembl"]
other_targets = targets.loc[lambda targets: targets['source'] == "other"]
gdrive_targets = targets.loc[lambda targets: targets['source'] == "gdrive"]

def get_sequence(wildcards, type):
	"""
	Return path to download sequence file
	"""
	species_units = targets.loc[wildcards.species]
	if type == "ensembl":
		if	species_units["source"] == "ensembl":
			return species_units["file"]
	if type == "tsa":
		if	species_units["source"] == "tsa":
			return species_units["file"]
	if type == "other":
		if species_units["source"] == "other":
			return species_units["file"]
	if type == "gdrive":
		if species_units["source"] == "gdrive":
			return species_units["file"]

rule get_ensembl:
	output:
		"resources/sequences/ensembl/{species}.fasta.gz" # this renames the original files and keep only the species names as given in the species column of targets dataframe
	params:
		lambda wc: get_sequence(wc, type="ensembl")
	shell:
		"rsync -v rsync://ftp.ensembl.org/ensembl/{params} {output}"

rule wget_tsa:
	output:
		"resources/sequences/tsa/{species}.fasta.gz"
	params:
		lambda wc: get_sequence(wc, type="tsa")
	shell:
		"wget -O resources/sequences/tsa/{wildcards.species}.fasta.gz {params}"

rule wget_other:
	output:
		"resources/sequences/other/{species}/{filename}"
	params:
		lambda wc: get_sequence(wc, type="other")
	shell:
		"""
		wget {params} -O {output}
		"""

rule wget_gdrive:
	output:
		"resources/sequences/gdrive/{species}.pep.fasta"
	params:
		lambda wc: get_sequence(wc, type="gdrive")
	shell:
		"wget -O {output} {params}" #filename and path is embedded in the gdrive command listed in the file column of target tsv
