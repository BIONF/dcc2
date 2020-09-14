# dcc2

DCCv2 is a tool for compiling core set data for [h1s](httsp://github.com/BIONF/HaMStR) using predicted orthologs from OMA, both OMA-browser and OMA-standalone. Outputs of this tool are 3 (optional 4) folders required for a HaMStR run, including (1) core_orthologs (comprises of OMA orthologous group - OG, or OMA pairs - OP. Each OG/OP has its own directory, where a multiple fasta file and a corresponding profile HMM can be found), (2) genome_dir (contains gene sets of taxa, from which the orthologs are originated), (3) blast_dir (holds the blast databases of those gene sets within genome_dir), and an optional (4) weight_dir (contains feature architecure annotations of all gene sets).

For OMA-standalone, DCCv2 requires the output orthoXML file from OMA, a taxon mapping file in tab-delimited format containing 3 columns (blabla) (blabla) (blabla). Protein set of included taxa can be either given as a folder, or automatically downloaded from OMA database.

To use DCCv2 with OMA-browser, some data from OMA need to be downloaded and processed in advance using the function Get OMA data under the DCCv2 menu of phylosophy. For OMA-browser, DCCv2 accept input as NCBI taxonomy IDs, taxon names or an OMA group ID. If only 2 taxa are given, one can decided if OMA groups or OMA pairs should be taken into account.
