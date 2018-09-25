import mygene
import pandas as pd
import re
import os
import numpy as np
import collections as cx
import os.path


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
    print(len(go_names)) # Includes ONLY signficant results
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

    if len(G) < 1:
        return 0, "NA"

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





def get_all(f1,f2,out, Dataset,Treatment,Species,Celltype,Selection,Windowsize, Coverage):
    df_final = pd.DataFrame(columns=["Dataset"], index=range(0,1))

    df_final["Dataset"] = Dataset
    df_final["Treatment"] = Treatment
    df_final["Species"] = Species
    df_final["Celltype"] = Celltype
    df_final["Selection"] = Selection
    df_final["Windowsize"] = Windowsize
    df_final["Coverage"] = Coverage
    
    
    def col_to_string(df):
        if len(df) > 0:
            s = ",".join(df["gene_id_name"].tolist())
            return s
        else:
            return "NA"



    df = pd.read_csv(f1,sep="\t")



    df_dog, df_pog = get_DOGS(df, "DOG")
    df_dogup, df_dogdo = get_sig_runion_csv(df_dog,0.05, log2FoldChange=0)
    df_pogup, df_pogdo = get_sig_runion_csv(df_pog,0.05, log2FoldChange=0)

    dfa = pd.read_csv(f2,sep="\t")

    dfa_dog, dfa_pog = get_DOGS(dfa, "ADOG")
    dfa_dogup, dfa_dogdo = get_sig_runion_csv(dfa_dog,0.05, log2FoldChange=0)
    dfa_pogup, dfa_pogdo = get_sig_runion_csv(dfa_pog,0.05, log2FoldChange=0)




    df_final["DOGs total"] = len(df_dog)                
    df_final["POGs total"] = len(df_pog)                
    df_final["ADOGs total"] = len(dfa_dog)              
    df_final["APOGs total"] = len(dfa_pog)              
    df_final["DOGs total names"] = col_to_string(df_dog)
    df_final["POGs total names"] = col_to_string(df_pog)
    df_final["ADOGs total names"] = col_to_string(dfa_dog)
    df_final["APOGs total names"] = col_to_string(dfa_pog)
    
    
    dfsig_dog  = get_sig_only_csv(df_dog, 0.05)
    dfsig_pog  = get_sig_only_csv(df_pog, 0.05)
    dfsiga_dog = get_sig_only_csv(dfa_dog, 0.05)
    dfsiga_pog = get_sig_only_csv(dfa_pog, 0.05)
    
    df_final["DOGs significant total"] = len(dfsig_dog)
    df_final["POGs significant total"] = len(dfsig_pog)
    df_final["ADOGs significant total"] = len(dfsiga_dog)
    df_final["APOGs significant total"] = len(dfsiga_pog)
    df_final["DOGs significant total names"] = col_to_string(dfsig_dog)
    df_final["POGs significant total names"] = col_to_string(dfsig_pog)
    df_final["ADOGs significant total names"] = col_to_string(dfsiga_dog)
    df_final["APOGs significant total names"] = col_to_string(dfsiga_pog)
    
    

    
    
    df_final["ADOGs up"] = len(dfa_dogup)               
    df_final["APOGs up"] = len(dfa_pogup)               
    df_final["DOGs up"] = len(df_dogup)                 
    df_final["POGs up"] = len(df_pogup)                 
    df_final["ADOGs down"] = len(dfa_dogdo)             
    df_final["APOGs down"] = len(dfa_pogdo)             
    df_final["DOGs down"] = len(df_dogdo)               
    df_final["POGs down"] = len(df_pogdo)               
    
    
    
    
    
    
    df_final["DOGs up names"] = col_to_string(df_dogup)
    df_final["POGs up names"] = col_to_string(df_pogup)
    df_final["ADOGs up names"] =    col_to_string(dfa_dogup)
    df_final["APOGs up names"] =    col_to_string(dfa_pogup)
    df_final["DOGs down names"] = col_to_string(df_dogdo)
    df_final["POGs down names"] = col_to_string(df_pogdo)
    df_final["ADOGs down names"] = col_to_string(dfa_dogdo)
    df_final["APOGs down names"] = col_to_string(dfa_pogdo)




    df_dog____total_genes, df_dog____line =  genes_to_GO(df_dog,colname="gene_id_name", outname=out + "All_DOG"    ,species=Species)
    df_pog____total_genes, df_pog____line =  genes_to_GO(df_pog,colname="gene_id_name", outname=out + "All_POG"    ,species=Species)
    dfa_dog___total_genes, dfa_dog___line =  genes_to_GO(dfa_dog,colname="gene_id_name", outname=out +"All_ADOG"   ,species=Species)
    dfa_pog___total_genes, dfa_pog___line =  genes_to_GO(dfa_pog,colname="gene_id_name", outname=out +"All_APOG"   ,species=Species)



    dfsig_dog____total_genes, dfsig_dog____line =  genes_to_GO(dfsig_dog,colname="gene_id_name", outname=out   +"sig_DOG"   ,species=Species)
    dfsig_pog____total_genes, dfsig_pog____line =  genes_to_GO(dfsig_pog,colname="gene_id_name", outname=out   +"sig_POG"   ,species=Species)
    dfsiga_dog___total_genes, dfsiga_dog___line =  genes_to_GO(dfsiga_dog,colname="gene_id_name", outname=out  +"sig_ADOG"  ,species=Species)
    dfsiga_pog___total_genes, dfsiga_pog___line =  genes_to_GO(dfsiga_pog,colname="gene_id_name", outname=out  +"sig_APOG"  ,species=Species)
    


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


    df_final["DOGs total Amount GO translation"]       =   df_dog____total_genes
    df_final["POGs total Amount GO translation"]       =   df_pog____total_genes
    df_final["ADOGs total Amount GO translation"]       =  dfa_dog___total_genes
    df_final["APOGs total Amount GO translation"]       =  dfa_pog___total_genes


    df_final["DOGs significant total Amount GO translation"]       =   dfsig_dog____total_genes
    df_final["POGs significant total Amount GO translation"]       =   dfsig_pog____total_genes
    df_final["ADOGs significant total Amount GO translation"]       =  dfsiga_dog___total_genes
    df_final["APOGs significant total Amount GO translation"]       =  dfsiga_pog___total_genes

    df_final["DOGs up Amount GO translation"]      =   len(list_match(dfsig_dog____line, df_dogup))
    df_final["POGs up Amount GO translation"]      =   len(list_match(dfsig_pog____line, df_pogup))
    df_final["ADOGs up Amount GO translation"]     =   len(list_match(dfsiga_dog___line, dfa_dogup))
    df_final["APOGs up Amount GO translation"]     =   len(list_match(dfsiga_pog___line, dfa_pogup))
    df_final["DOGs down Amount GO translation"]    =   len(list_match(dfsig_dog____line, df_dogdo))
    df_final["POGs down Amount GO translation"]    =   len(list_match(dfsig_pog____line, df_pogdo))
    df_final["ADOGs down Amount GO translation"]   =   len(list_match(dfsiga_dog___line, dfa_dogdo))
    df_final["APOGs down Amount GO translation"]   =   len(list_match(dfsiga_pog___line, dfa_pogdo))



    df_final["DOGs total GO translation"]   =      df_dog____line
    df_final["POGs total GO translation"]   =      df_pog____line
    df_final["ADOGs total GO translation"]  =      dfa_dog___line
    df_final["APOGs total GO translation"]  =      dfa_pog___line

    df_final["DOGs up GO translation"]                   =       ",".join(list_match(dfsig_dog____line, df_dogup)  )
    df_final["POGs up GO translation"]                   =       ",".join(list_match(dfsig_pog____line, df_pogup)  )
    df_final["ADOGs up GO translation"]                  =       ",".join(list_match(dfsiga_dog___line, dfa_dogup) )
    df_final["APOGs up GO translation"]                  =       ",".join(list_match(dfsiga_pog___line, dfa_pogup) )
    df_final["DOGs down GO translation"]                 =       ",".join(list_match(dfsig_dog____line, df_dogdo)  )
    df_final["POGs down GO translation"]                 =       ",".join(list_match(dfsig_pog____line, df_pogdo)  )
    df_final["ADOGs down GO translation"]                =       ",".join(list_match(dfsiga_dog___line, dfa_dogdo) )
    df_final["APOGs down GO translation"]                =       ",".join(list_match(dfsiga_pog___line, dfa_pogdo) )

    return df_final









#Worms


HSPATH = "../data/Dogcatcher_Out/HS/N2_vs_HS_2ndtime/FINAL_OUT/ALL"
OKPATH = "../data/Dogcatcher_Out/OK/N2_vs_OK_2ndtime/FINAL_OUT/ALL"



f1 =  HSPATH +"/DOG/combined_DOG_noOperon.csv"
f2 =  HSPATH +"/ADOG/combined_ADOG_noOperon.csv"
df_HS = get_all(f1,f2,"GO_terms_DOGS/melnick2018/HS",Dataset="melnick2018",Treatment="heat shock 36 °C 3 hr",Species="worm",Celltype="whole worm", Selection="J2 IP", Windowsize=100, Coverage=90)
# #
f1 =  OKPATH +"/DOG/combined_DOG_noOperon.csv"
f2 =  OKPATH +"/ADOG/combined_ADOG_noOperon.csv"
df_OK = get_all(f1,f2,"GO_terms_DOGS/melnick2018/OK",Dataset="melnick2018",Treatment="tdp-1 deletion(ok803)",Species="worm",Celltype="whole worm", Selection="J2 IP", Windowsize=100, Coverage=90)





df_final = pd.concat([df_HS, df_OK])
df_final.to_csv("Meta_comparison_DOGS.tsv",sep="\t", index=None)
