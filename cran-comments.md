


checking Rd cross-references ... WARNING
Unknown package 'herit' in Rd xrefs

checking for missing documentation entries ... WARNING
Undocumented data sets:
  'tr_cap_phenos_met'
All user-level objects in a package should have documentation entries.
See chapter 'Writing R documentation files' in the 'Writing R
Extensions' manual.



Documented arguments not in \usage in documentation object 'herit':
  'ms_exp'

Undocumented arguments in documentation object 'herit_boot'
  'ms_exp'

Functions with \usage entries need to have the appropriate \alias
entries, and all their arguments documented.
The \usage entries must correspond to syntactically valid R code.
See chapter 'Writing R documentation files' in the 'Writing R
Extensions' manual.

checking DESCRIPTION meta-information ... NOTE
Package listed in more than one of Depends, Imports, Suggests, Enhances:
  'purrr'
A package should be listed in only one of these fields.

checking top-level files ... NOTE
Non-standard files/directories found at top level:
  'README.Rmd' 'other_data' 'scripts'

checking R code for possible problems ... NOTE
dist_env: no visible global function definition for 'combn'
dist_env: no visible global function definition for 'as.dist'
herit.lm: no visible global function definition for 'anova'
herit.lm: no visible binding for global variable 'term'
herit.lm: no visible binding for global variable 'meansq'
herit.lmerMod: no visible binding for global variable 'grp'
herit.lmerMod: no visible binding for global variable 'vcov'
herit_boot.lm: no visible global function definition for 'model.frame'
herit_boot.lm: no visible global function definition for 'simulate'
... 8 lines ...
herit_boot.lmerMod: no visible global function definition for 'sd'
herit_boot.lmerMod: no visible global function definition for
  'quantile'
Undefined global functions or variables:
  anova as.dist combn grp map_dbl meansq model.frame quantile sd
  simulate term update vcov
Consider adding
  importFrom("stats", "anova", "as.dist", "model.frame", "quantile",
             "sd", "simulate", "update", "vcov")
  importFrom("utils", "combn")
to your NAMESPACE file.

  
  
  
herit.lm: no visible binding for global variable ‘term’ (C:\Users\Jeff\Documents\Programming\pbr/R/herit.R:116-117)
herit.lm: no visible binding for global variable ‘meansq’ (C:\Users\Jeff\Documents\Programming\pbr/R/herit.R:116-117)
herit.lmerMod: no visible binding for global variable ‘grp’ (C:\Users\Jeff\Documents\Programming\pbr/R/herit.R:165-166)
herit_boot.lm: no visible global function definition for ‘map_dbl’ (C:\Users\Jeff\Documents\Programming\pbr/R/herit_boot.R:111)
herit_boot.lm: local variable ‘mf’ assigned but may not be used (C:\Users\Jeff\Documents\Programming\pbr/R/herit_boot.R:89)
herit_boot.lmerMod: local variable ‘mf’ assigned but may not be used (C:\Users\Jeff\Documents\Programming\pbr/R/herit_boot.R:142)
