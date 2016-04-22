# This function is hidden function of GOSemSim package
gene2GO <- function(gene, organism="human", ont="MF", dropCodes="IEA") {
  library(GOSemSim)
  gene <- as.character(gene)
  loadGOMap(organism)
  gomap <- get("gomap", envir=GOSemSimEnv)
  go <- gomap[[gene]]
  
  if (all(is.na(go)))
    return (NA)
  
  ## go.df <- ldply(go, function(i) c(GOID=i$GOID, Evidence=i$Evidence, Ontology=i$Ontology))
  ## go.df <- go.df[ !go.df$Evidence %in% dropCodes, ] ## compatible to work with NA and NULL
  ## goid <- go.df[go.df$Ontology == ont, "GOID"]
  goid <- sapply(go, function(i) i$GOID)
  evidence <- sapply(go, function(i) i$Evidence)
  ontology <- sapply(go, function(i) i$Ontology)
  
  idx <- ! evidence %in% dropCodes
  goid <- goid[idx] ## drop dropCodes Evidence
  ontology <- ontology[idx]
  goid <- goid[ontology == ont]
  
  if (length(goid) == 0)
    return (NA)
  
  goid <- as.character(unique(goid))
  return (goid)
}
