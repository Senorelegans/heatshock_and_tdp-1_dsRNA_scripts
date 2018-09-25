"""Data downloaded from NCBI Gene data converted into Python data."""

# Put this link in then modify it for your species of interest

#From goatools.......
# http://www.ncbi.nlm.nih.gov/gene/?term=genetype+protein+coding%5BProperties%5D+AND+%229606%22%5BTaxonomy+ID%5D+AND+alive%5Bproperty%5D
# or follow link to "Search for Human Protein-Coding Genes" found here:
# https://github.com/dvklopfenstein/biocode/blob/master/biodownloads/README.md


import collections as cx
import pandas as pd

df = pd.read_csv("genes_worm_ncbi_download.txt", sep="\t",index_col=False)

names = " ".join(list(df))

NtData = cx.namedtuple('NtData', names)

# NtData = cx.namedtuple('NtData', 'tax_id Org_name GeneID CurrentID Status Symbol Aliases description other_designations map_location chromosome genomic_nucleotide_accession_version start_position_on_the_genomic_accession end_position_on_the_genomic_accession orientation exon_count OMIM no_hdr0')

GeneID2nt = {}

for x in range(0,len(df)):
    id = df["GeneID"].iloc[x]
    # print(id)
    GeneID2nt[id] = NtData(
        tax_id=df["tax_id"].iloc[x],
        Org_name=df["Org_name"].iloc[x],
        GeneID=df["GeneID"].iloc[x],
        CurrentID=df["CurrentID"].iloc[x],
        Status=df["Status"].iloc[x],
        Symbol=df["Symbol"].iloc[x],
        Aliases=df["Aliases"].iloc[x],
        description=df["description"].iloc[x],
        other_designations=df["other_designations"].iloc[x],
        map_location=df["map_location"].iloc[x],
        chromosome=df["chromosome"].iloc[x],
        genomic_nucleotide_accession_version=df["genomic_nucleotide_accession_version"].iloc[x],
        start_position_on_the_genomic_accession=df["start_position_on_the_genomic_accession"].iloc[x],
        end_position_on_the_genomic_accession=df["end_position_on_the_genomic_accession"].iloc[x],
        orientation=df["orientation"].iloc[x],
        exon_count=df["exon_count"].iloc[x],
        OMIM=df["OMIM"].iloc[x],
    )


# GENEID2NT = {
#   1 : NtData(tax_id=9606, Org_name='Homo sapiens', GeneID=1, CurrentID=0, Status='live', Symbol='A1BG', Aliases='A1B, ABG, GAB, HYST2477', description='alpha-1-B glycoprotein', other_designations='HEL-S-163pA|epididymis secretory sperm binding protein Li 163pA', map_location='19q13.4', chromosome='19', genomic_nucleotide_accession_version='NC_000019.10', start_position_on_the_genomic_accession=58346806, end_position_on_the_genomic_accession=58353499, orientation='minus', exon_count=8, OMIM=138670, no_hdr0=''),
#   2 : NtData(tax_id=9606, Org_name='Homo sapiens', GeneID=2, CurrentID=0, Status='live', Symbol='A2M', Aliases='A2MD, CPAMD5, FWP007, S863-7', description='alpha-2-macroglobulin', other_designations='C3 and PZP-like alpha-2-macroglobulin domain-containing protein 5|alpha-2-M', map_location='12p13.31', chromosome='12', genomic_nucleotide_accession_version='NC_000012.12', start_position_on_the_genomic_accession=9067708, end_position_on_the_genomic_accession=9115962, orientation='minus', exon_count=36, OMIM=103950, no_hdr0='')}
