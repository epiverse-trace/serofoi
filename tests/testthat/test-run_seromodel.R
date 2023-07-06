set.seed(123)
data("serodata")
serodata <- prepare_serodata(serodata)
foi_constant <- run_seromodel(serodata, foi_model = "constant")
