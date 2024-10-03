import mygene

def map_human_ensembl_to_gene_symbols(ensembl):
    try:
        mg_sess = map_human_ensembl_to_gene_symbols.session
    except AttributeError:
        mg_sess = mygene.MyGeneInfo()
        map_human_ensembl_to_gene_symbols.session = mg_sess
    mapped = mg_sess.querymany(
        ensembl,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        verbose=False,
        returnall=True,
    )

    # build a mapping of the genes between symbols and ensembl
    gene_map = {}
    for g in mapped["out"]:
        q = g["query"]
        try:
            sym = g["symbol"]
        except KeyError:
            sym = None
        gene_map[q] = sym
    return gene_map