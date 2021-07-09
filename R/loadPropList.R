

loadProp = function(p_prop){
  
  # load prop
  d_in = read.csv(p_prop, sep = "\t", row.names = 1)
  drops = c("Group", "Aff", "SMILES")
  d_prop = d_in[ , !(names(d_in) %in% drops)]
  
  # for E2up
  d_prop$E2up[which(is.na(d_prop$E2up))] = "NT"
  d_prop$E2up[which(d_prop$E2up == "ns effect")] = "NEG"
  d_prop$E2up[which(d_prop$E2up != "NEG" & d_prop$E2up != "NT")] = "POS"
  
  #for P4up
  d_prop$P4up[which(is.na(d_prop$P4up))] = "NT"
  d_prop$P4up[which(d_prop$P4up == "ns effect")] = "NEG"
  d_prop$P4up[which(d_prop$P4up != "NEG" & d_prop$P4up != "NT")] = "POS"
  
  #for ER
  d_prop$ER[which(d_prop$ER == "")] = "NT"
  d_prop$ER[which(d_prop$ER == "tested")] = "NEG"
  d_prop$ER[which(d_prop$ER == "agonist" | d_prop$ER == "antagonist")] = "POS"
  
  #genotox
  d_prop$genotox[which(d_prop$genotox == "not tested" | d_prop$genotox == "predicted genotoxic" | d_prop$genotox == "" | is.na(d_prop$genotox))] = "NT"
  d_prop$genotox[which(d_prop$genotox == "not genotoxic")] = "NEG"
  d_prop$genotox[which(d_prop$genotox != "NEG" & d_prop$genotox != "NT")] = "POS"
  
  #forMC
  d_prop$MC[which(d_prop$MC == "" | is.na(d_prop$MC))] = "NT"
  d_prop$MC[which(d_prop$MC == "0")] = "NEG"
  d_prop$MC[which(d_prop$MC == "1")] = "POS"
  
  #H295R
  d_prop$H295R[which(d_prop$H295R == "" | is.na(d_prop$H295R))] = "NT"
  d_prop$H295R[which(d_prop$H295R == "0")] = "NEG"
  d_prop$H295R[which(d_prop$H295R == "1")] = "POS"
  
  return(d_prop)
  
}
