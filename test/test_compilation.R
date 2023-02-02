remotes::install_github("TRACE-LAC/serofoi", ref = "dev")

library(serofoi)

data_test <- prepare_data(mydata)

model_0_object <- run_model(
  model_data = data_test,
  model_name = "constant_foi_bi",
  n_iters = 1000
)
model_0_plot <- plot_model(model_0_object, size_text = 6)

plot_seroprev(model_0_object, size_text = 15)
plot_foi(model_0_object, size_text = 15)
plot_rhats(model_0_object, size_text = 15)
summary_model <- extract_summary_model(model_0_object)
