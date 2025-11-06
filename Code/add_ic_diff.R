# This function takes in the output dataframe from run_ic_folder() and adds the 
# differences in AIC and BIC between different orders 
# Arguments:
  # ic_df: output dataframe from run_ic_folder
# Output: 
  # Dataframe containing differences in AIC and BIC between orders

add_ic_diff <- function(ic_df) {
  num_cols <- grep("^(AIC|BIC)[0-5]$", names(ic_df), value = TRUE)
  ic_df[num_cols] <- lapply(ic_df[num_cols], function(x) as.numeric(as.character(x)))
  
  # Add AIC and BIC differences
  ic_df$dAIC_2_minus_1 <- ic_df$AIC2 - ic_df$AIC1
  ic_df$dBIC_2_minus_1 <- ic_df$BIC2 - ic_df$BIC1
  ic_df$dAIC_1_minus_0 <- ic_df$AIC1 - ic_df$AIC0
  ic_df$dBIC_1_minus_0 <- ic_df$BIC1 - ic_df$BIC0

  return(ic_df)
}