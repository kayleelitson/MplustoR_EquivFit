#title: "Functions for Equivalence Testing Simulation Study: R Import of Mplus Output Files"
#author of R code: "BLINDED FOR PEER REVIEW" (some code taken directly from Marcoulides & Yuan functions, 2017)
#other authors on paper: "BLINDED FOR PEER REVIEW"
#code last updated: "5/9/2022"

# The code here it to pull in Mplus input files to then estimate the T-size RMSEA & CFI statistics described
# by the article "Confirm Structural Equation Models by Equivalence Testing with Adjusted Fit Indices"
# by Yuan, Chan, Marcoulides and Bentler. We Use MplusAutomation by Michael Hallquist & Wiley (2018)
# to pull in output files for different types of modeling approaches. We removed cases with missing values
# on the fit statistics necessary for computation and recomputed necessary models by noting the number of missing
# per condition. We auto remove duplicates and have functions for removing files with missing fit from analysis.



##################### BEGIN FUNCTION FOR NOTING WHICH PIECES OF OUTPUT ARE NOT MISSING AND CAN BE USED #####################
whichOutComplete = function(outFile){
  outFile <- outFile[unlist(lapply(outFile, length) != 0)]
  data0 <- NULL;
  data1 <- NULL;
  for (i in 1:length(outFile)){
    data0$filename[i] <- outFile[[i]]$summaries$Filename %>% as.character() %>% tolower();
  }
  data0$filename <- unique(data0$filename)
  
  #identification specific to these output files!
  data0 <- as.data.frame(data0);
  data0 <- separate(data0, filename, 
                    into = c("model", "varsHL", "misspec", "missing", "skew", "n_pop", "iteration"),
                    sep = "_", remove = FALSE);
  data0 <- separate(data0, iteration, into = c("iteration", NA), sep = "[.]", remove = TRUE);
  col_names <- names(data0);
  data0[ , col_names] <- lapply(data0[ , col_names] , factor);
  #end of specific to these output files.
  
  data0$ID = 1:nrow(data0) %>% as.integer;
  data0$miss <- NULL;
  data0$miss <- 1:nrow(data0) %>% as.integer;
  data0$warn <- NULL;
  data0$warn <- 1:nrow(data0) %>% as.integer;
  data0 <- as.data.frame(data0);
  
  ###Exclude missing files from analysis. Print missing files & count of missing files by ID vars in separate output
  for (i in 1:nrow(data0)){
    if("ChiSqM_Value" %in% names(outFile[[i]]$summaries)
       & "ChiSqM_DF" %in% names(outFile[[i]]$summaries)
       & "ChiSqBaseline_Value" %in% names(outFile[[i]]$summaries)
       & "ChiSqBaseline_DF" %in% names(outFile[[i]]$summaries)
       & "NDependentVars" %in% names(outFile[[i]]$summaries)
       & "ChiSqM_PValue" %in% names(outFile[[i]]$summaries)){
      data0$miss[i] = 0 
    } else {
      data0$miss[i] = 1
    }
  }
  
  ###Note which outputs have warnings
  for (i in 1:nrow(data0)){
    if(length(outFile[[i]][["warnings"]]) > 0) {
      data0$warn[i] = 1 
    } else {
      data0$warn[i] = 0
    }
  }
  
  data0 <- as.data.frame(data0);
  data0$miss <- as.integer(data0$miss);
  data0$warn <- as.integer(data0$warn); #NEW
  data1 <- data0[ which(data0$miss == 0), ]
  #for (i in 1:nrow(data0)){outFileNoMissing <- outFile[ which(outFile[[i]] == data1$ID)]}
  return(data1);
}


##################### BEGIN FUNCTION FOR SAVING ONLY COMPLETE OUTPUT THAT CAN BE USED #####################
saveComplete = function(outFile){
  outFile <- outFile[unlist(lapply(outFile, length) != 0)]
  data0 <- NULL;
  data1 <- NULL;
  outFileNoMissing <- NULL;
  for (i in 1:length(outFile)){
    outFile[[i]]$ID <- i %>% as.integer()
    data0$filename[i] <- outFile[[i]]$summaries$Filename %>% as.character() %>% tolower()
    data0$ID[i] <- outFile[[i]]$ID
  }
  data0$filename <- unique(data0$filename)
  
  #identification specific to these output files!
  data0 <- as.data.frame(data0);
  data0 <- separate(data0, filename, 
                    into = c("model", "varsHL", "misspec", "missing", "skew", "n_pop", "iteration"),
                    sep = "_", remove = FALSE);
  data0 <- separate(data0, iteration, into = c("iteration", NA), sep = "[.]", remove = TRUE);
  col_names <- names(data0);
  data0[ , col_names] <- lapply(data0[ , col_names] , factor);
  #end of specific to these output files.
  
  data0$miss <- NULL;
  data0$miss <- 1:nrow(data0) %>% as.integer;
  data0 <- as.data.frame(data0);
  
  ###Exclude missing files from analysis. Print missing files & count of missing files by ID vars in separate output
  for (i in 1:nrow(data0)){
    if("ChiSqM_Value" %in% names(outFile[[i]]$summaries)
       & "ChiSqM_DF" %in% names(outFile[[i]]$summaries)
       & "ChiSqBaseline_Value" %in% names(outFile[[i]]$summaries)
       & "ChiSqBaseline_DF" %in% names(outFile[[i]]$summaries)
       & "NDependentVars" %in% names(outFile[[i]]$summaries)
       & "ChiSqM_PValue" %in% names(outFile[[i]]$summaries)){
      data0$miss[i] = 0 
    } else {
      data0$miss[i] = 1
    }
  }
  
   ###Note which outputs have warnings
  for (i in 1:nrow(data0)){
    if(length(outFile[[i]][["warnings"]]) > 0) {
      data0$warn[i] = 1 
    } else {
      data0$warn[i] = 0
    }
  }
  
  data0 <- as.data.frame(data0);
  data0$miss <- as.integer(data0$miss);
  data0$warn <- as.integer(data0$warn);

  data0$n_samp <- NA
  data0$chisq <- NA
  data0$degf <- NA
  data0$chipval <- NA
  data0$chisq_null <- NA
  data0$degf_null <- NA
  data0$p_obsvars <- NA

  for (i in 1:nrow(data0)){
    if(data0$miss[i] == 0 & data0$warn == 0){
      data0$n_samp[i] <- outFile[[i]]$summaries[["Observations"]]
      data0$chisq[i] <- outFile[[i]]$summaries[["ChiSqM_Value"]]
      data0$degf[i] <- outFile[[i]]$summaries[["ChiSqM_DF"]]
      data0$chipval[i] <- outFile[[i]]$summaries[["ChiSqM_PValue"]]
      data0$chisq_null[i] <- outFile[[i]]$summaries[["ChiSqBaseline_Value"]]
      data0$degf_null[i] <- outFile[[i]]$summaries[["ChiSqBaseline_DF"]]
      data0$p_obsvars[i] <- as.integer(outFile[[i]]$summaries[["NDependentVars"]]) + as.integer(outFile[[i]]$summaries[["NIndependentVars"]])
    } else {
      data0$n_samp[i] <- NA
      data0$chisq[i] <- NA
      data0$degf[i] <- NA
      data0$chipval[i] <- NA
      data0$chisq_null[i] <- NA
      data0$degf_null[i] <- NA
      data0$p_obsvars[i] <- NA
    }
  }
  
  data1 <- data0[ which(data0$miss == 0 & data0$warn == 0), ]
  data1 <- as.data.frame(data1);

  return(data1)
}



##################### BEGIN FUNCTION FOR NCP_CHI2, NECESSARY FOR EQUIVALENCE TEST #####################
# The formula for ncp_chi2 is from Venables 1975 for obtaining the noncentrality 
# of a non-central chi-square distribution, resulting in ncp_chi2 value in the function below.
ncp_chi2=function(alpha, T_ml, df){
  z=qnorm(1-alpha);
  z2=z*z; z3=z2*z; z4=z3*z; z5=z4*z;
  sig2=2*(2*T_ml-df+2);
  sig2 = ifelse(sig2 < 0, 0, sig2); #this line was added to constrain negative sig2 values to 0,
                                    #otherwise an error occurred in the sqrt fx in the line below.
  sig=sqrt(sig2); sig3=sig*sig2; sig4=sig2*sig2;sig5=sig4*sig;
  sig6=sig2*sig4;
  
  delta=T_ml-df+2+sig*
    (
      z+(z2-1)/sig-z/sig2 + 2*(df-1)*(z2-1)/(3*sig3)
      +( -(df-1)*(4*z3-z)/6+(df-2)*z/2 )/sig4
      +4*(df-1)*(3*z4+2*z2-11)/(15*sig5)
      +(
        -(df-1)*(96*z5+164*z3-767*z)/90-4*(df-1)*(df-2)*(2*z3-5*z)/9
        +(df-2)*z/2
      )/sig6
    );
  delta=max(delta, 0, na.rm = TRUE);
  return(delta)
}



##################### BEGIN FUNCTION FOR EQUIVALENCE TESTING #####################

modelFit_equivalence_test = function(df_outFile=dataframe_of_completeMplusOutput, alpha=alpha){
  ###Create & Define Variable Names
  data1 <- df_outFile;
  alpha <- .05;
  # filename = name of the output file
  # date = date
  # model = which model?
  # vars_hl = High or low number of variables in input
  # misspec = degree of misspecification
  # missing = degree of missingness
  # skew = degree of skewness
  # N_pop = sample size of the population model
  # iteration = iteration (1 - 5000)
  # N_samp = sample size of the sample (after accounting for missing)
  # chisq = Chi-square model fit value
  # degf = chi-square degrees of freedom value
  # chipval = chi-square p value for the model
  # chisq_null = Chi-square null value
  # degf_null = chi-square degrees of freedom null value
  # p_obsvars = total number of observed variables in the model
  
  
  # Creating Variables to input to the Marcoulides & Yuan Equations
  # Input and Calculating RMSEA_t and CFI_t: needed inputs are the observed statistic T_ml, 
  # its degrees of freedom (df), sample size (N). For estimating T-size CFI, additional inputs 
  # are the observed statistic at the independence model T_mli, and the number of observed variables (p);
  # Updated 7-28-2021 to include na.rm = TRUE for all max() functions
  for (i in 1:nrow(df_outFile)){
    data1$degf_i[i] = data1$p_obsvars[i]*(data1$p_obsvars[i]+1)/2 - data1$p_obsvars[i];
    data1$delta_c[i] = max(0, data1$chisq[i] - data1$degf[i], na.rm = TRUE);
    data1$delta_i[i] = data1$chisq_null[i] - data1$degf_i[i];
    data1$delta_tR[i] = ncp_chi2(alpha, data1$chisq[i], data1$degf[i]); #delta_t for RMSEA
    data1$delta_tC[i] = ncp_chi2(alpha/2, data1$chisq[i], data1$degf[i]); #delta_t for CFI
    data1$delta_it[i] = ncp_chi2(1 - alpha/2, data1$chisq_null[i], data1$degf_i[i]); #delta_it
    data1$RMSEA_c[i] = sqrt(data1$delta_c[i] / (data1$degf[i]*(data1$n_samp[i] - 1)) );
    data1$RMSEA_t[i] = sqrt(data1$delta_tR[i] / (data1$degf[i]*(data1$n_samp[i] - 1)) );
    data1$CFI_c[i] = 1 - data1$delta_c[i] / max(data1$delta_c[i], data1$delta_i[i], 0.0000000001, na.rm = TRUE);
    data1$CFI_t[i] = 1 - max(data1$delta_tC[i], 0, na.rm = TRUE) / max(data1$delta_tC[i], data1$delta_it[i], 0.0000000001, na.rm = TRUE);
    #computing cut values
    data1$n[i] = data1$n_samp[i] - 1; 
    data1$degf_i[i] = data1$p_obsvars[i]*(data1$p_obsvars[i] - 1)/2;
    data1$CFI_e99[i] = 1-exp( 4.67603-.50827*log(data1$degf[i])+.87087*(data1$degf[i]^(1/5))-.59613*((data1$degf_i[i])^(1/5))-1.89602*log(data1$n[i])
                              + .10190*((log(data1$n[i]))^2)+ .03729*log(data1$degf[i])*log(data1$n[i]) );
    data1$CFI_e95[i] = 1-exp( 4.12132-.46285*log(data1$degf[i])+.52478*(data1$degf[i]^(1/5))-.31832*((data1$degf_i[i])^(1/5))-1.74422*log(data1$n[i])
                              +.13042*((log(data1$n[i]))^2)-.02360*(data1$n[i]^(1/2))+.04215*log(data1$degf[i])*log(data1$n[i]) );
    data1$CFI_e92[i] = 1-exp( 6.31234-.41762*log(data1$degf[i])+.01554*((log(data1$degf[i]))^2)-.00563*((log(data1$degf_i[i]))^2)-1.30229*log(data1$n[i])
                              +.19999*((log(data1$n[i]))^2)-2.17429*(data1$n[i]^(1/5))+.05342*log(data1$degf[i])*log(data1$n[i])-.01520*log(data1$degf_i[i])*log(data1$n[i]) );
    data1$CFI_e90[i] = 1-exp( 5.96633-.40425*log(data1$degf[i])+.01384*((log(data1$degf[i]))^2)-.00411*((log(data1$degf_i[i]))^2)-1.20242*log(data1$n[i])
                              +.18763*((log(data1$n[i]))^2)-2.06704*(data1$n[i]^(1/5))+.05245*log(data1$degf[i])*log(data1$n[i])-.01533*log(data1$degf_i[i])*log(data1$n[i]) );
    data1$RMSEA_e01[i] = exp( 1.34863-.51999*log(data1$degf[i])+.01925*log(data1$degf[i])*log(data1$degf[i])-.59811*log(data1$n[i])+.00902*sqrt(data1$n[i])
                              +.01796*log(data1$degf[i])*log(data1$n[i]) );
    data1$RMSEA_e05[i] = exp(2.06034-.62974*log(data1$degf[i])+.02512*log(data1$degf[i])*log(data1$degf[i])-.98388*log(data1$n[i])+.05442*log(data1$n[i])*log(data1$n[i])
                             -.00005188*data1$n[i]+.05260*log(data1$degf[i])*log(data1$n[i]) );
    data1$RMSEA_e08[i] = exp( 2.84129-.54809*log(data1$degf[i])+.02296*log(data1$degf[i])*log(data1$degf[i])-.76005*log(data1$n[i])+.10229*log(data1$n[i])*log(data1$n[i])
                              -1.11167*(data1$n[i]^.2)+.04845*log(data1$degf[i])*log(data1$n[i]) );
    data1$RMSEA_e10[i] = exp( 2.36352-.49440*log(data1$degf[i])+.02131*log(data1$degf[i])*log(data1$degf[i])-.64445*log(data1$n[i])+.09043*log(data1$n[i])*log(data1$n[i])
                              -1.01634*(data1$n[i]^.2)+.04422*log(data1$degf[i])*log(data1$n[i]) );
    
    #Output for conventional and t-size CFI and RMSEA values & their cut values#
    data1$RMSEA_conv[i] = data1$RMSEA_c[i];
    data1$RMSEA_t[i] = data1$RMSEA_t[i];
    #data1$RMSEA_cutoff01[i] = round(data1$RMSEA_e01[i], digits = 3) #Lower Bound of Excellent Fit
    data1$RMSEA_cutoff05[i] = round(data1$RMSEA_e05[i], digits = 3) #Lower Bound of Good Fit. Lowest Acceptable Model Fit
    #data1$RMSEA_cutoff08[i] = round(data1$RMSEA_e08[i], digits = 3) #Lower Bound of Fair Fit
    #data1$RMSEA_cutoff10[i] = round(data1$RMSEA_e10[i], digits = 3) #Lower Bound of Mediocre Fit. Below this value is Poor Fit.
    data1$CFI_conv[i] = data1$CFI_c[i];
    data1$CFI_t[i] = data1$CFI_t[i];
    #data1$CFI_cutoff99[i] = round(data1$CFI_e99[i], digits = 3) #Lower Bound of Excellent Fit
    data1$CFI_cutoff95[i] = round(data1$CFI_e95[i], digits = 3) #Lower Bound of Good Fit. Lowest Acceptable Model Fit
    #data1$CFI_cutoff92[i] = round(data1$CFI_e92[i], digits = 3) #Lower Bound of Fair Fit
    #data1$CFI_cutoff90[i] = round(data1$CFI_e90[i], digits = 3) #Lower Bound of Mediocre Fit. Below this value is Poor Fit.
  }
  
  data1 = as.data.frame(data1);
  
  for (i in 1:nrow(df_outFile)){
    #goodfit = 1 for good and excellent fit (above cutoff95 value for CFI, below cutoff05 value for RMSEA, and non-significant (> .05) for chi-square p value)
    data1$RMSEA_t_goodfit[i][data1$RMSEA_t[i] <= data1$RMSEA_cutoff05[i]] <- 1
    data1$RMSEA_t_goodfit[i][data1$RMSEA_t[i] > data1$RMSEA_cutoff05[i]] <- 0
    data1$RMSEA_t_goodfit[i][is.na(data1$RMSEA_t[i])] <- NA
    
    data1$CFI_t_goodfit[i][data1$CFI_t[i] >= data1$CFI_cutoff95[i]] <- 1
    data1$CFI_t_goodfit[i][data1$CFI_t[i] < data1$CFI_cutoff95[i]] <- 0
    data1$CFI_t_goodfit[i][is.na(data1$CFI_t[i])] <- NA
    
    data1$chisq_goodfit[i][data1$chipval[i] > 0.0500] <- 1
    data1$chisq_goodfit[i][data1$chipval[i] <= 0.0500] <- 0
    data1$chisq_goodfit[i][is.na(data1$chipval[i])] <- NA
  }
  
  #Create variables of for fit according to correctly and incorrectly specified models
    #fit_match_cond = 1 for good and excellent fit and = 0 for fair and worse fit when the model is CORRECTLY specified. 
    #fit_match_cond = 0 for good and excellent fit and = 1 for fair and worse fit when the model is INCORRECTLY specified.
    data1 <- mutate(data1, RMSEA_t_fit_match_cond = ifelse(grepl("correct", misspec),
                                                           data1$RMSEA_t_goodfit,
                                                           abs(data1$RMSEA_t_goodfit - 1)),
                    CFI_t_fit_match_cond = ifelse(grepl("correct", misspec),
                                                  data1$CFI_t_goodfit,
                                                  abs(data1$CFI_t_goodfit - 1)),
                    chisq_fit_match_cond = ifelse(grepl("correct", misspec),
                                                  data1$chisq_goodfit,
                                                  abs(data1$chisq_goodfit - 1)))
  
  data1 = as.data.frame(data1);

  data2 <- subset(data1, select=c("filename", "model", "varsHL", "misspec", "missing", "skew", "n_pop", "iteration",
                  "n_samp", "chisq", "degf", "chipval", "chisq_null", "degf_null", "p_obsvars",
                  "RMSEA_conv", "RMSEA_t", "RMSEA_cutoff05", "CFI_conv","CFI_t", "CFI_cutoff95",
                  "chisq_goodfit", "RMSEA_t_goodfit", "CFI_t_goodfit",
                  "chisq_fit_match_cond", "RMSEA_t_fit_match_cond", "CFI_t_fit_match_cond"
                  ))
  return(data2); #change to data1 for all variables, data2 for subset
}

#Function to Restrict data to 5000 cases per condition (for cases where the number of replications exceeded 5000)
restrict5000 <- function(dataframe1){
  dataframe1$count <- NULL
  dataframe1$group <- NULL
  dataframe1 <- dataframe1 %>% mutate(group = group_indices(., .dots= 
                                                               c("model", "varsHL", 
                                                                 "misspec", "missing", 
                                                                 "skew", "n_pop")))
  dataframe1$count <- ave(dataframe1$group, dataframe1$group, FUN=seq_along)
  dataframe1$count <- as.numeric(dataframe1$count)
  dataframe1 <- dataframe1[which(dataframe1$count <= 5000),]
  return(dataframe1)
}


