import mygene
import pandas as pd
import re
import os
import numpy as np
import collections as cx
from goatools.godag_plot import plot_goid2goobj
# Get http://geneontology.org/ontology/go-basic.obo# Get h
# from goatools.base import download_go_basic_obo
# obo_fname = download_go_basic_obo()
#
# from goatools.base import download_ncbi_associations
# gene2go = download_ncbi_associations()

from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GeneID2nt as GeneID2nt_human
# from goatools.test_data.genes_NCBI_6239_ProteinCoding import GeneID2nt as GeneID2nt_worm
from goatools.go_enrichment import GOEnrichmentStudy

from genes_NCBI_6239_marko_ALL import GeneID2nt as GeneID2nt_worm



def genes_to_GO(df,colname, outname, species):
    GeneID2nt_temp = {}
    taxid = 0

    if species == "human":
        taxid = 9606
        GeneID2nt_temp = GeneID2nt_mus
    if species == "mouse":
        taxid = 10090
        GeneID2nt_temp = GeneID2nt_human
    if species == "worm":
        taxid = 6239
        GeneID2nt_temp = GeneID2nt_worm

    if len(df) < 1:
        return 0, "NA"
    mg = mygene.MyGeneInfo()
    genelist = df[colname].tolist()
    G = []
    for gene in genelist:
        if ("ENS") in gene:
            if ("T") in gene:
                tcl = [m.start() for m in re.finditer(r"\.",gene)]  #Tab count list
                gene = gene[:tcl[0]]
        G.append(gene)

    df = pd.DataFrame(mg.querymany(G, scopes='ensembl.gene,ensembl.transcript,symbol,reporter,accession,retired,hgnc,entrezgene,name,wormbase', fields='entrezgene', species=taxid, as_dataframe=True))

    # Remove ensembl or refseq without id's
    print("Amount queryed: ", len(df))
    df = df[df["_id"].str.contains("_") == False ]
    df = df[df["_id"].str.contains("ENS") == False ]
    df = df[df["_id"].str.contains("WB") == False ]
    df["symbol"] = df.index
    df["id"] = df["_id"]
    df = df[["symbol","id"]]
    print("Amount kept with id: ", len(df))
    if len(df) < 1:
        return 0, "NA"

    ######################## GOA TOOLS ########################
    obodag = GODag("go-basic.obo")
    geneid2gos_species = read_ncbi_gene2go("gene2go", taxids=[taxid])
    print(len(geneid2gos_species) ,  "annotated genes")

    if species == "human":
        goeaobj = GOEnrichmentStudy(
            GeneID2nt_human.keys(), # List of species protein-coding genes
            geneid2gos_species, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method
    if species == "mouse":
        goeaobj = GOEnrichmentStudy(
            GeneID2nt_mus.keys(), # List of species protein-coding genes
            geneid2gos_species, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method
    if species == "worm":
        goeaobj = GOEnrichmentStudy(
                GeneID2nt_worm.keys(), # List of species protein-coding genes
                geneid2gos_species, # geneid/GO associations
                obodag, # Ontologies
                propagate_counts = False,
                alpha = 0.05, # default significance cut-off
                methods = ['fdr_bh']) # defult multipletest correction method

    # Data will be stored in this variable# Data
    geneid2symbol = {}
    for i in range(len(df)):
        geneid2symbol[int(df["id"].iloc[i])] = df["symbol"].iloc[i]

    # 'p_' means "pvalue". 'fdr_bh' is the multipletest method we are currently using.# 'p_'
    geneids_study = geneid2symbol.keys()
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

    ######################################
    #Change study_items back to gene names
    goeaobj.wr_xlsx(outname + ".xlsx", goea_results_sig)
    def studyitemToGene(genes):
        if type(genes) == str:
            genes = genes.split(",")
            genes = [geneid2symbol[int(g)] for g in genes]
            if species == "worm":
                df = pd.DataFrame(mg.querymany(genes, scopes='ensembl.gene,ensembl.transcript,symbol,reporter,accession,retired,hgnc,entrezgene,name,wormbase', fields='symbol', species=taxid, as_dataframe=True))
                genes = df["symbol"].tolist()
                print(genes)
            genes = ",".join(genes)
            print("returned genes are ", genes)
            return genes
        else:
            return ""

    if os.path.isfile(outname + ".xlsx"):
        df = pd.read_excel(outname + ".xlsx", index=None)
        df["gene_id_names"] = df.apply(lambda row: studyitemToGene(row["study_items"]), axis=1)
        del df["study_items"]
        print(df["gene_id_names"])
        df.to_excel(outname + ".xlsx")
    ######################################

    go_names = [r.name for r in goea_results_sig]
    # print(len(go_names)) # Includes ONLY signficant results
    word2cnt = cx.Counter([word for name in go_names for word in name.split()])
    # Print 10 most common words found in significant GO term names
    print(word2cnt.most_common(10))
    freq_seen = ['RNA', 'translation', 'mitochond', 'ribosomal', 'ribosome']
    freq_seen = ['translation','ribosomal', 'ribosome']

    # Collect significant GOs for words in freq_seen (unordered)
    word2siggos = cx.defaultdict(set)
    # Loop through manually curated words of interest
    for word in freq_seen:
        # Check each significant GOEA result for the word of interest
        for rec in goea_results_sig:
            if word in rec.name:
                word2siggos[word].add(rec.GO)
    # Sort word2gos to have the same order as words in freq_seen
    word2siggos = cx.OrderedDict([(w, word2siggos[w]) for w in freq_seen])

    goid2goobj_all  = {nt.GO:nt.goterm for nt in goea_results_all}
    go2res = {nt.GO:nt for nt in goea_results_all}

    fout = outname + "_GO_word_genes.txt"
    total_genes = 0
    G = []

    for word, gos in word2siggos.items():
        # Sort first by BP, MF, CC. Sort second by GO id.
        gos = sorted(gos, key=lambda go: [go2res[go].NS, go])
        genes = set()
        for go in gos:
            genes |= go2res[go].study_items
        genes = sorted([geneid2symbol[g] for g in genes])
        G = G + genes
        # print("\n{WD}: {N} study genes, {M} GOs\n".format(WD=word, N=len(genes), M=len(gos)))
        # print("{WD} GOs: {GOs}\n".format(WD=word, GOs=", ".join(gos)))
        for i, go in enumerate(gos):
            res = go2res[go]
            print("{I}) {NS} {GO} {NAME} ({N} genes)\n".format(
                I=i, NS=res.NS, GO=go, NAME=res.name, N=res.study_count))


        N = 10 # Number of genes per line
        mult = [genes[i:i+N] for i in range(0, len(genes), N)]

    G = list((set(G))) # Remove duplicates
    total_genes = len(G)

    G = ",".join(G)
    print("total genes are", total_genes)
    print("G is", G)

    return total_genes, G


#### Process dogs from datasets

def get_DOGS(df, DOG):
    if DOG == "DOG":
        df_DOG = df[df["TYPE"]=="DOG"]
        df_POG = df[df["TYPE"]=="POG"]
    if DOG == "ADOG":
        df_DOG = df[df["TYPE"]=="ADOG"]
        df_POG = df[df["TYPE"]=="APOG"]
    return df_DOG, df_POG


def get_sig_runion_csv(df, padj, log2FoldChange):
    """This function will take in a DSeq2 normalized matrix and make a gtf of all non-sig genes"""
    df = df.copy()
    df["significant"] = np.where(df["padj"] < padj, True, False)
    df_final = df[df.significant].copy()  #Keep all True values. (significant)
    del df_final["significant"]
    df_final = df_final.copy()
    #Treatment vs Control
    df_final["T_vs_C_up"] = np.where(df_final["log2FoldChange"] > log2FoldChange, True, False)
    df_T_vs_C_up = df_final[df_final.T_vs_C_up]  #Keep all True values. (up)
    del df_T_vs_C_up["T_vs_C_up"]
    df_final["T_vs_C_down"] = np.where(df_final["log2FoldChange"] < -1*log2FoldChange, True, False)
    df_T_vs_C_down = df_final[df_final.T_vs_C_down]  #Keep all True values. (down)
    del df_T_vs_C_down["T_vs_C_up"]
    del df_T_vs_C_down["T_vs_C_down"]
    return df_T_vs_C_up, df_T_vs_C_down


def get_sig_only_csv(df, padj):
    """This function will take in a DSeq2 normalized matrix and make a gtf of all non-sig genes"""
    df = df.copy()
    df["significant"] = np.where(df["padj"] < padj, True, False)
    df_final = df[df.significant].copy()  #Keep all True values. (significant)
    del df_final["significant"]
    df_final = df_final.copy()
    return df_final



def col_to_string(df, col):
    if len(df) > 0:
        s = ",".join(df[col].tolist())
        return s
    else:
        return "NA"


def get_all(fs,fa,foverlap, f2, f3, out, Dataset,Treatment,Species,Celltype,Selection,Windowsize, Coverage):
    df_final = pd.DataFrame(columns=["Dataset"], index=range(0,1))

    df_final["Dataset"] = Dataset
    df_final["Treatment"] = Treatment
    df_final["Species"] = Species
    df_final["Celltype"] = Celltype
    df_final["Selection"] = Selection
    df_final["Windowsize"] = Windowsize
    df_final["Coverage"] = Coverage

    dfsense = pd.read_csv(fs,sep=",")
    dfantis = pd.read_csv(fa,sep=",")
    dfsense.rename(columns={"Unnamed: 0" : "gene_id_name"}, inplace=True)
    dfantis.rename(columns={"Unnamed: 0" : "gene_id_name"}, inplace=True)
    dfsenseup, dfsensedo = get_sig_runion_csv(dfsense,0.05, log2FoldChange=0)
    dfantisup, dfantisdo = get_sig_runion_csv(dfantis,0.05, log2FoldChange=0)


    dfsigsense = get_sig_only_csv(dfsense,0.05)
    dfsigantis = get_sig_only_csv(dfantis,0.05)
    
    

    df_final["Sense significant total amount"] = len(dfsigsense)
    df_final["Antisense significant total amount"] = len(dfsigantis)

    df_final["Sense significant total"] = col_to_string(dfsigsense, "gene_id_name")
    df_final["Antisense significant total"] = col_to_string(dfsigantis, "gene_id_name")




    df_final["Sense genes up"]          =  len(dfsenseup)
    df_final["Sense genes down"]        =  len(dfsensedo)

    df_final["Sense genes up names"]          =  col_to_string(dfsenseup, "gene_id_name")
    df_final["Sense genes down names"]        =  col_to_string(dfsensedo, "gene_id_name")


    dfsenseup_total_genes, dfsenseup_line =  genes_to_GO(dfsenseup,colname="gene_id_name", outname=out +"sense_up"   ,species=Species)
    df_final["Sense genes up Amount GO translation"] =      dfsenseup_total_genes
    df_final["Sense genes up GO translation"] =  dfsenseup_line
    dfsensedo_total_genes, dfsensedo_line =  genes_to_GO(dfsensedo,colname="gene_id_name", outname=out +"sense_down"   ,species=Species)
    df_final["Sense down Amount GO term translation"] =      dfsensedo_total_genes
    df_final["Sense down GO translation"] =  dfsensedo_line


    df_final["Antisense genes amount up"]      =  len(dfantisup)
    df_final["Antisense genes amount down"]    =  len(dfantisdo)

    df_final["Antisense genes up"]      =  col_to_string(dfantisup, "gene_id_name")
    df_final["Antisense genes down"]    =  col_to_string(dfantisdo, "gene_id_name")



    dfantisup_total_genes, dfantisup_line =  genes_to_GO(dfantisup,colname="gene_id_name", outname=out +"antisense_up"   ,species=Species)
    df_final["Antisense genes up amount translation"] =      dfantisup_total_genes
    df_final["Antisense genes up translation"] =  dfantisup_line
    dfantisdo_total_genes, dfantisdo_line =  genes_to_GO(dfantisdo,colname="gene_id_name", outname=out +"antisense_down"   ,species=Species)
    df_final["Antisense genes down amount translation"] =      dfantisdo_total_genes
    df_final["Antisense genes down translation"] =  dfantisdo_line


    #Get all ADOGs
    dfa = pd.read_csv(f2,sep="\t")
    dfa_dog, dfa_pog = get_DOGS(dfa, "ADOG")
    dfadogup, dfadogdo = get_sig_runion_csv(dfa,0.05, log2FoldChange=0)
    
    #Get DOGS
    dfdog = pd.read_csv(f3,sep="\t")
    df_dog, df_pog = get_DOGS(dfdog, "DOG")
    dfdog.rename(columns={"Unnamed: 0" : "gene_id_name"}, inplace=True)
    dfdogup, dfdogdo = get_sig_runion_csv(dfdog,0.05, log2FoldChange=0)
    

    #Get overlapped genes on opposite strand and rename with gene_id_name so you can add to ADOGs
    df_o = pd.read_csv(foverlap,sep="\t")
    df_o = df_o[df_o["same_or_opposite_strand"]=="opposite"]
    df_overlapped = df_o[["overlapped_gene"]]
    df_overlapped = df_overlapped.drop_duplicates(subset="overlapped_gene")

    df_final["Antisense opposite strand overlapped genes total"]  = len(df_overlapped)
    df_final["Antisense opposite strand overlapped genes total names"]  = col_to_string(df_overlapped, "overlapped_gene")



    #Get sig dog and adogs with sig overlapped antisense genes
    #Get antisense genes with sig dogs    
    dfdogup_w_overlap =  df_o[df_o["gene_id_name"].isin(dfdogup["gene_id_name"])]
    dfdogdo_w_overlap =  df_o[df_o["gene_id_name"].isin(dfdogdo["gene_id_name"])]
    dfantisup_overlapped =  dfdogup_w_overlap[dfdogup_w_overlap["overlapped_gene"].isin(dfantisup["gene_id_name"])]
    dfantisdo_overlapped =  dfdogdo_w_overlap[dfdogdo_w_overlap["overlapped_gene"].isin(dfantisdo["gene_id_name"])]
    dfantisup_overlapped["gene_id_name"] = dfantisup_overlapped["overlapped_gene"]  #make concat easier
    dfantisdo_overlapped["gene_id_name"] = dfantisdo_overlapped["overlapped_gene"]  #make concat easier
    
    
    #Get antisense genes with sig adogs
    dfantisup_adog =  dfadogup[dfadogup["gene_id_name"].isin(dfantisup["gene_id_name"])] 
    dfantisdo_adog =  dfadogdo[dfadogdo["gene_id_name"].isin(dfantisdo["gene_id_name"])] 
    
    
    
    
    
    dfup_OA = pd.concat([dfantisup_overlapped,dfantisup_adog])    # Combine overlapped and ADOGs
    dfup_OA = dfup_OA.drop_duplicates(subset="gene_id_name")
    dfdo_OA = pd.concat([dfantisdo_overlapped,dfantisdo_adog])    # Combine overlapped and ADOGs
    dfdo_OA = dfdo_OA.drop_duplicates(subset="gene_id_name")
    

    df_final["Antisense opposite strand overlapped genes amount up"]      =  len(dfup_OA)
    df_final["Antisense opposite strand overlapped genes amount down"]    =  len(dfdo_OA)

    df_final["Antisense opposite strand overlapped genes up"]      =  col_to_string(dfup_OA, "gene_id_name")
    df_final["Antisense opposite strand overlapped genes down"]    =  col_to_string(dfdo_OA, "gene_id_name")

    dummy_var2, dummy_var=  genes_to_GO(dfup_OA,colname="gene_id_name", outname=out +"Antisense_opposite_strand_overlapped_genes_up"   ,species=Species)
    dummy_var2, dummy_var=  genes_to_GO(dfdo_OA,colname="gene_id_name", outname=out +"Antisense_opposite_strand_overlapped_genes_down"   ,species=Species)


    def list_match(s, df):
            print("List match ****************************************")
            list1 = s.split(",")
            list2 = df["gene_id_name"].tolist()
            G = []
            for gene in list2:
                if ("ENS") in gene:
                    if ("T") in gene:
                        tcl = [m.start() for m in re.finditer(r"\.",gene)]  #Tab count list
                        gene = gene[:tcl[0]]
                G.append(gene)

            print(list1)
            print(G)
            listmerged = list(set(list1).intersection(G))

            print(listmerged)
            return listmerged


    #dfantisup_line,dfantisup_overlapped
    df_final["Antisense opposite strand overlapped genes up amount translation"] =   len(list_match(dfantisup_line,dfup_OA))
    df_final["Antisense opposite strand overlapped genes down amount translation"] = len(list_match(dfantisdo_line,dfdo_OA))
    df_final["Antisense opposite strand overlapped genes up translation"] =        ",".join(list_match(dfantisup_line,dfup_OA))
    df_final["Antisense opposite strand overlapped genes down translation"] =      ",".join(list_match(dfantisdo_line,dfdo_OA))
    #


    return df_final


#Worm

HSPATH = "../data/Dogcatcher_Out/HS/initial_Rsubread"
HSDOGPATH = "../data/Dogcatcher_Out/HS/N2_vs_HS_2ndtime/FINAL_OUT/ALL"

fs = HSPATH + "/DESeq2_sense.csv"
fa = HSPATH + "/DESeq2_antisense.csv"
foverlap = HSDOGPATH + "/DOG/biotypes/combined_DOG_with_biotypes_UNPACKED.csv"
f2 = HSDOGPATH + "/ADOG/combined_ADOG_noOperon.csv"
f3 = HSDOGPATH + "/DOG/combined_DOG_noOperon.csv"
df_HS = get_all(fs,fa,foverlap, f2, f3, "GO_terms_genes/melnick2018/HS",Dataset="melnick2018",Treatment="heat shock 36 °C 3 hr",Species="worm",Celltype="whole worm", Selection="J2 IP", Windowsize=100, Coverage=90)




OKPATH = "../data/Dogcatcher_Out/OK/initial_Rsubread"


OKDOGPATH = "../data/Dogcatcher_Out/OK/N2_vs_OK_2ndtime/FINAL_OUT/ALL"
fs = OKPATH + "/DESeq2_sense.csv"
fa = OKPATH + "/DESeq2_antisense.csv"
foverlap = OKDOGPATH + "/DOG/biotypes/combined_DOG_with_biotypes_UNPACKED.csv"
f2 = OKDOGPATH + "/ADOG/combined_ADOG_noOperon.csv"
f3 = OKDOGPATH + "/DOG/combined_DOG_noOperon.csv"
df_OK = get_all(fs,fa,foverlap, f2, f3, "GO_terms_genes/melnick2018/OK",Dataset="melnick2018",Treatment="tdp-1 deletion(ok803)",Species="worm",Celltype="whole worm", Selection="J2 IP", Windowsize=100, Coverage=90)






df_final = pd.concat([df_HS, df_OK])
# df_final = pd.concat([df_HS])



df_final.to_csv("Meta_comparison_antisense_overlapped.tsv",sep="\t", index=None)

