
## Debugging codes

source("r/function/MakeADFun_windows_debug.R")
library(TMB)
TMB::compile("models/3_reallocation_discrete_toy_model.cpp")
dyn.load("models/3_reallocation_discrete_toy_model")
MakeADFun_windows_debug(cpp_name = "models/3_reallocation_discrete_toy_model",
                        data=Data, parameters=Params,  random=Random)
TMB::gdbsource("models/3_reallocation_discrete_toy_model.R",interactive = T) ## Non-interactive
dyn.unload( dynlib("models/3_reallocation_discrete_toy_model" ) )
