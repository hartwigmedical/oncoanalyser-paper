
rules <- 
list(
  "events" = list(
     "tso" = list(
        "miss" = c("fusion", "hla", "viral", "homozygous disruption", "hrd"),
        "capture" = c("tml high", "tmb high","msi high")),
      "amp" = list(
         "miss" = c("fusion", "hla", "viral", "homozygous disruption", "hrd", "tml high", "tmb high"),
         "capture" = c("msi high"))  
  ), 
  "cn" = list("min" = .5, "max" = 6)
)

apply_coverage_rules <- function( event, rules = rules, panel = "tso", val ){
    if( tolower(event) %in% rules$events[[panel]]$miss){ FALSE}
    else if (tolower(event) %in% rules$events[[panel]]$capture) { TRUE }
    else { grepl("TRUE", val)}
}
