# MultiCens

MultiCens is platform from where we can predict hormone-responsive genes in a target tissue by providing the inputs such as Source and target tissue names, the mediator hormone, and the list of genes from source and the target tissues. The tool will take the list of genes from the source tissue and will rank the hormone-responsive genes from the target tissue gene list by their responsiveness for the selected hormone. And we can also measure the Centrality scores such as local centrality, global centrality, and query-set centrality.

Local Centrality
A node in a layer (i.e. Tissue) can affect other nodes (i.e. genes) in the same layer as well as different layers. Local centrality captures the within-layer effect of a node.

Global Centrality
Local centrality considers the effect of only within-layer connections, Global centrality captures the remaining effect. The global centrality of a node is a measure of its influence on all nodes irrespective of their layers.

Query Set Centrality
Query-set centrality captures the effect of a node on a query-set of nodes present in any specific layer in the multilayer network.


Using Tool v2.0, We can measure the centrality scores for genes such as local centrality, global centrality, and query-set centrality by using gene expression data from GTex dataset or by providing your own gene expression data for Tissues.
                       
