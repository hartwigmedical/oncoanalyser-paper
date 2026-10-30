level_map <- list(
  "FDA Approved"                  = "Rest FDA",
  "Drug Repurposing"              = "Rest Repurpose",
  "Experimental" = "Rest Experimental",
  "Other"                         = "Rest Other",
  "None"                          = "None")

name_map <- list(
  "FDA Approved"                  = "FDA",
  "Drug Repurposing"              = "Off",
  "Experimental" = "Exp",
  "Other"                         = "Oth",
  "None"                          = "None")

highest_level_event <- function(a, b, c, d) {
  if      (a > 0) "FDA"
  else if (b > 0) "Off"
  else if (c > 0) "Exp"
  else if (d > 0) "Oth"
  else            "None"
}

panel_maps <- list(
  "wgs"     = "WGS",
  "tso_cov" = "Comprehensive\n521 gene panel",
  "amp_cov" = "Targeted    \n50 gene panel")
