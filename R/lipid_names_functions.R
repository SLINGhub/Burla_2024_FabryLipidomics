### Load R libraries

library(dplyr)
library(stringr)
library(tidyr)





###################################################################################################################
# Clean Lipid Names
#  - removes "C" in front of carbon numbers
#  TODO: The rest we DONT NEED AND WANT....
###################################################################################################################


lipid_name_cleaner <- function(lipid_name){

  name_temp <- str_replace(lipid_name, "\\/C(?=\\d)","\\/")
  name_temp <- str_replace(name_temp, "[:space:]C"," ")

  # name_temp <- str_replace(name_temp, "dhCer ","Cer d18:0/")
  # name_temp <- str_replace(name_temp, "DG 15:0 15:0 (IS)","DG 30:0 [-15:0] (IS)")
  # name_temp <- str_replace(name_temp, "DG 38:6_NL 20:5","DG 38:6 -(20:5)")
  # name_temp <- str_replace(name_temp, "DG 40:6_NL 22:5","DG 40:6 -(22:5)")
  # name_temp <- str_replace(name_temp, "DG 40:7_NL 20:5","DG 40:7 -(22:6)")
  # name_temp <- str_replace(name_temp, "DG 40:6_NL 22:6","DG 40:6 -(22:6)")
  # name_temp <- str_replace(name_temp, "DG 40:8_NL 22:6","DG 40:8 -(22:6)")
  # name_temp <- str_replace(name_temp, "17:0 17:0","34:0")
  # name_temp <- str_replace(name_temp, "17:0/17:0","34:0")
  # name_temp <- str_replace(name_temp, "13:0 13:0","26:0")
  # name_temp <- str_replace(name_temp, "-d"," d")
  # name_temp <- str_replace(name_temp, "\\(104\\)"," [104]")
  # name_temp <- str_replace(name_temp, "\\(161\\)"," [161]")
  # TODO:  Forgot what this does. adding ] to DG TG transitions with FA?
  #if(str_detect(name_temp, fixed("DG|TG"))){
  #  toreplace <- str_extract(name_temp, "\\-\\(.*?\\)")
  #  replaceby <- str_replace(toreplace,"\\)","]")
  #  name_temp <- str_replace(name_temp, fixed(toreplace), replaceby)
  #}
  name_temp <- str_replace(name_temp, "\\-\\(","[-")
  name_temp <- str_replace(name_temp, "\\/C","/")
  name_temp <- str_replace(name_temp, "\\((?=[0-9, d, m, t, O, P])"," ")
  name_temp <- str_replace(name_temp, "(?<=[0-9])\\)+","")
  name_temp <- str_replace(name_temp, "  "," ")

  return (name_temp)
}


###################################################################################################################
# Get and Add Transition and Compound-Specific Info
###################################################################################################################

get_transition_info <- function(transitions, add_sum_composition, add_chain_info, add_sum_composition_name = TRUE){

  temp_transitions <- transitions %>% add_lipid_transition_classnames()

  if(add_sum_composition | add_chain_info){
    chain_info <- transitions %>%  get_lipid_sumcomposition(add_sum_composition, add_chain_info, add_sum_composition_name)
    temp_transitions <- temp_transitions %>% left_join(chain_info, by=c("Compound"="Compound"))

    if(add_sum_composition_name){
      temp_transitions <- temp_transitions |>
        mutate(sum_composition_name = paste0(str_trim(lipidClass), " ", chain_length_total, ":", chain_doublebonds_total))
    }


  }


  #dat_transitions <- transitions %>% left_join(temp_transitions, by=c("Compound"="Compound"))
  return(temp_transitions)
}




###################################################################################################################
# Determine Lipid and Transition Class Name from Compound name
###################################################################################################################
# Retrieves "lipid class" from the compound name and adds it a factor column. Lipid class is defined as group of lipids sharing same head group and modications but with different chain lengths and saturations.
# lipidClassBase are for example all ceramides (including deoxy etc) and PC (including PC-O etc)
# Different transitions of the same lipid class are indicatec in transition class

get_CompoundName <- function(transition_name){
  compound_name <- str_trim(str_replace(transition_name, "(?<!TG) \\[.*?\\]",""))
  compound_name <- str_trim(str_replace(compound_name, "  "," "))
  compound_name <- str_trim(str_replace(compound_name, "\\+",""))
  compound_name <- str_trim(str_replace(compound_name, "(?<!Cer) \\[.*?\\]",""))
  compound_name <- str_trim(str_replace(compound_name, "  "," "))


  return(compound_name)
}

add_lipid_transition_classnames <- function(datLong){
  #cat("Retrieving lipid class/transition names...", fill = FALSE)
  datLong_temp <- datLong %>%  mutate(lipidClassBase = (str_trim(str_extract(Compound, "[A-z0-9]+[[:blank:]]*"))))
  #datLong_temp <- datLong_temp %>%  mutate(lipidClass = (str_trim(str_extract(Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m])"))))

  # add a "-", except for between name and sphingoid base
  #datLong_temp <- datLong_temp %>% mutate(lipidClass = str_replace(lipidClass, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"))
  datLong_temp <- datLong_temp %>% mutate(lipidClass = str_squish(ifelse(str_detect(Compound, "\\;"),str_c(str_extract(Compound, "[A-z0-9]+[[:blank:]]"), " ", str_extract(Compound, "(?<=\\;).{2}")), str_extract(Compound, "[A-z0-9]+[[:blank:]]"))))
  #datLong_temp <- datLong_temp %>% mutate(lipidClassSL = str_extract(Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9\\:]{4})"))
  datLong_temp <- datLong_temp %>% mutate(lipidClassSL = str_extract(Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[0-9\\:\\/\\;O]{7})"))
  # Add  transition name  (defined by flanking "[]")
  datLong_temp <- datLong_temp %>%  mutate(transitionName = str_extract(Compound, "(?<=\\[).*?(?=\\])"))
  datLong_temp <- datLong_temp %>%  mutate(transitionName = ifelse(is.na(transitionName),"",transitionName))

  # Add  transition class (adds measured product class in square brackets after the compound class name). In case of NL indicated with [-...] (e.g.  TG 48:2 [-16:1]), [NL] is indicated (TG [NL]). Maybe useful for normalization
  datLong_temp <- datLong_temp %>%  mutate(transitionClass =  ifelse(grepl("\\[\\-|\\[NL", Compound),"M>M-NL",str_replace(paste("", str_trim(str_extract(Compound, "\\[[a-zA-Z0-9\\-\\: ]*\\]")))," NA", transitionName)))

  datLong_temp <- datLong_temp %>%  mutate(CompoundName = get_CompoundName(Compound))
  #datLong_temp <- datLong_temp %>%  mutate(CompoundName =  str_trim(str_replace(Compound, "\\[.*?\\]","")))
  #datLong_temp <- datLong_temp %>%  mutate(CompoundName =  str_trim(str_replace(CompoundName, "  "," ")))

  #datLong_temp <- datLong_temp %>% mutate_at(.vars = c("lipidClass","lipidClassBase", "transitionClass", "transitionName" , "CompoundName"), as.factor)

  # Convert the transition names for each lipidClass to a number (e.g. for Cer d18:1 "M>SphB","M>SphB-H2O"  becomes 1,2:  for for Cer m18:1 "M-H2O>SphB", "M-H2O>SphB" also becomes 1, 2)
  datLong_temp <- datLong_temp %>%
    group_by(CompoundName) %>%
    do(
      mutate(.data, transition_group = match(.$transitionName , levels(as.factor(datLong_temp[datLong_temp$lipidClassSL == unique(.$lipidClassSL),]$transitionName))))
    ) %>%
    mutate(transition_group = ifelse(is.na(transition_group), 1, transition_group)) %>%
    mutate(CompoundName = if_else(str_detect(Compound, "^TG"), Compound, CompoundName)) %>%
    relocate(CompoundName, lipidClassBase, lipidClass, lipidClassSL, transitionClass, transitionName, transition_group, .after = Compound)

  datLong_temp <- datLong_temp |> mutate(lipidClassSL = if_else(str_detect(lipidClassSL, "S1P"), "S1P", lipidClassSL))
  datLong_temp <- datLong_temp |> mutate(lipidClassSL = if_else(str_detect(CompoundName, "LysoGb3"), "LysoGb3", lipidClassSL))

  datLong_temp <- datLong_temp |>
    mutate(lipidClassBase = if_else(str_detect(lipidClass, "\\-O|\\-P"), lipidClass, lipidClassBase))

  datLong_temp <- datLong_temp |>
    mutate(lipidClass = if_else(str_detect(CompoundName, "LysoGb3"), "LysoGb3 O2", lipidClass))

  datLong_temp <- datLong_temp |>
    mutate(lipidClassBase =
             case_when(
               str_detect(CompoundName, "Hex|GM|Gm") ~ "GSL",
               str_detect(CompoundName, "HexCer d18\\:0|Hex1Cer d18\\:0|Hex2Cer d18\\:0") ~ "dhGSL",
               str_detect(CompoundName, "Cer d18\\:0|Cer d16\\:0|Cer d17\\:0|Cer d19\\:0|Cer d20\\:0") ~ "dhCer",
               str_detect(CompoundName, "Cer d") ~ "Cer",
               str_detect(CompoundName, "Cer m") ~ "deoxyCer",
               str_detect(CompoundName, "Ga2") ~ "GSL",
               str_detect(CompoundName, "Cer m") ~ "deoxyCer",
               str_detect(CompoundName, "S1P|Sph") ~ "SPB",
               str_detect(CompoundName, "SM [1-2]") ~ "SM",
               str_detect(CompoundName, "SM [3-9]") ~ "SM_sum",
               str_detect(CompoundName, "LysoGb3") ~ "LGSL",
               str_detect(CompoundName, "LHex2Cer") ~ "LGSL",
               TRUE ~ lipidClassBase
             ))

  return(datLong_temp %>% ungroup())

}


###################################################################################################################
# Gets  FA and LCB chains from lipid names and calculates sum compositions (total Chain length and double bonds)
###################################################################################################################
# INFO: Takes max.3 chains (TAG). Does not work e.g. for cardiolipins etc
# TODO: HANDLE SIDE CHAINS !!!
# TODO: Robustness and validating results (is number etc)
get_lipid_sumcomposition <- function(d_compounds, add_sum_composition, add_chain_info, add_sum_composition_name = FALSE){

  d <- d_compounds %>% dplyr::select(Compound) %>% distinct()
  d_temp = d %>% mutate(chains = str_extract(Compound, "(?<=\\s).*?(?=(?:\\s|$))"))
  d_temp = d_temp %>% mutate(chains = str_replace(chains, "O\\-|P\\-|d|t|m",""))
  d_temp <- d_temp %>% separate(chains, into = c("chain_1", "chain_2", "chain_3"),sep = "/|_", remove = TRUE,fill = "right")
  d_temp <- d_temp %>% separate(chain_1, into = c("chain_1_C", "chain_1_DB"),sep = ":", remove = TRUE,fill = "right",convert=TRUE)
  d_temp <- d_temp %>% separate(chain_2, into = c("chain_2_C", "chain_2_DB"),sep = ":", remove = TRUE,fill = "right",convert=TRUE)
  d_temp <- d_temp %>% separate(chain_3, into = c("chain_3_C", "chain_3_DB"),sep = ":", remove = TRUE,fill = "right",convert=TRUE)
  d_temp <- d_temp %>% mutate_at(.vars=vars(matches("chain")),
                                 .funs = ~str_extract(., "\\-*\\d+\\.*\\d*"))
  d_temp <- d_temp %>% mutate_at(.vars=vars(matches("chain")),
                                 .funs= ~as.numeric(.))
  d_temp <- d_temp %>%
    mutate(across(everything(),~replace(., is.na(.), 0))) %>%
    mutate(chain_length_total = dplyr::select(., contains("_C")) %>% rowSums(),
           chain_doublebonds_total = dplyr::select(., contains("_DB")) %>% rowSums())
  if (!add_chain_info) d_temp <- d_temp %>% dplyr::select(Compound, chain_length_total, chain_doublebonds_total)


  if (!add_sum_composition) d_temp <- d_temp %>% dplyr::select(Compound, -chain_length_total, -chain_doublebonds_total)


  return(d_temp %>% ungroup())
}

clean_lipidnames <- function(datLong){
  datLong_clean <- datLong %>%
    mutate(Compound= str_replace(Compound, "\\/C(?=\\d)","\\/")) %>%
    mutate(Compound= str_replace(Compound, "[:space:]C"," ")) %>%
    mutate(Compound= str_replace(Compound, "(?<=[A-z])_(?=[0-9])"," ")) %>%
    mutate(Compound= str_replace(Compound, "(?<=[A-z])_(?=d[0-9])"," "))%>%
    mutate(Compound= str_replace(Compound, "\\(IS\\)"," (IS)")) %>%
    mutate(Compound= str_replace(Compound, "d7"," d7")) %>%
    mutate(Compound= str_replace(Compound, "d9"," d9"))%>%
    mutate(Compound= str_replace(Compound, "PC-O ","PC O-")) %>%
    mutate(Compound= str_replace(Compound, "PC-P ","PC P-")) %>%
    mutate(Compound= str_replace(Compound, "PE-O ","PE O-")) %>%
    mutate(Compound= str_replace(Compound, "PE-P ","PE P-")) %>%
    mutate(Compound= str_replace(Compound, "SM d","SM ")) %>%
    mutate(Compound= str_replace(Compound, "Cer d","Cer ")) %>%
    mutate(Compound= str_replace(Compound, "SM 36-2","SM 36:2"))
  return(datLong_clean)
}

# Retreive lipid class from lipid name
add_lipidclass_names <- function(datLong){
  datLong_temp <- clean_lipidnames(datLong)
  datLong_temp <- datLong_temp %>%
    mutate(lipidClassBase = str_trim(str_extract(Compound, "[A-z0-9]+[[:blank:]]*")),
           lipidClass = str_trim(str_extract(Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9]{2}[:]{1}[0-9]{1})")))

  datLong_temp <- datLong_temp %>%
    mutate(lipidClass = str_replace(lipidClass, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"),
           isPUFA = ifelse(str_detect(Compound, "\\:0|\\:1|\\:2"),FALSE,TRUE))

  datLong_temp <- datLong_temp %>% mutate(lipidClass = ifelse(lipidClassBase=="Cer", "Cer", lipidClass))
  datLong_temp <- datLong_temp %>% mutate(lipidClass = ifelse(lipidClassBase=="Cer", "Cer", lipidClass))
  datLong_temp <- datLong_temp %>% mutate(lipidClass = ifelse(lipidClassBase=="HexCer", "HexCer", lipidClass))
  datLong_temp <- datLong_temp %>% mutate_at(.vars = c("lipidClass", "lipidClassBase"), as.factor)


  return(datLong_temp)
}



add_chains_sum <- function(datLong){

  get_chain_composition <- function(Compound){
    chain_description = str_extract(Compound, "(?<=\\s).*?(?=(?:\\s|$))")
    chain_description = str_replace(chain_description, "O\\-|P\\-|d|t|m","")
    chains <- str_split(chain_description,"/|_")
    chain_details = str_split(unlist(chains),":", simplify = TRUE)
    class(chain_details) <- "numeric"
    comp_sum <- colSums(chain_details)
    return(tibble(chain_length_total = comp_sum[1],
                  chain_doublebonds_total = comp_sum[2]))

  }

  datLong_temp <- datLong %>%  mutate(chain_info = map(Compound, get_chain_composition))  %>% unnest(chain_info)
  return(datLong_temp)
}

