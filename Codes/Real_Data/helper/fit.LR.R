###########################################################################
# Presence/Absence Test using Logistic Regression (Only returns p-values) #
###########################################################################

#######################
# Fit LR To A Dataset #
#######################

fit.LR <- function(features,
                   metadata,
                   libSize,
                   ID){

  #####################
  # Per-feature model #
  #####################

  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){

    ###############################
    # Extract features one by one #
    ###############################

    featuresVector <- features[, x]
    featuresVector<-ifelse(featuresVector!=0, 1, 0) # Convert to Presence/Absence

    #################################
    # Create per-feature input data #
    #################################

    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize, ID)
    formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))

    ##############################################
    # Automatic library size adjustment for GLMs #
    ##############################################

    if(length(unique(libSize)) > 1){ # To prevent offsetting with TSS-normalized data
      formula<-update(formula, . ~ . - offset(log(libSize)))
    }

    ########################
    # Random effects model #
    ########################

    if(!length(ID) == length(unique(ID))){
      formula<-update(formula, . ~ . + (1|ID))
      fit <- tryCatch({
        fit1 <- lme4::glmer(formula = formula,
                            data = dat_sub,
                            family = "binomial")

      }, error=function(err){
        fit1 <- try({lme4::glmer(formula = formula,
                                 data = dat_sub,
                                 family = "binomial")})
        return(fit1)
      })

      ###################################
      # Summarize Coefficient Estimates #
      ###################################

      if (class(fit) != "try-error"){
        para<-as.data.frame(summary(fit)$coefficients)[-1,-3]
      } else{
        print(paste("Fitting problem for feature", x, "returning NA"))
        para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=3))
      }
      colnames(para)<-c('coef_LR', 'stderr_LR', 'pval_LR')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
    }

    #######################
    # Fixed effects model #
    #######################

    else{
      fit <- tryCatch({
        fit1 <- glm(formula, data = dat_sub, family = 'binomial')
      }, error=function(err){
        fit1 <- try({glm(formula, data = dat_sub, family = 'binomial')})
        return(fit1)
      })

      ###################################
      # Summarize Coefficient Estimates #
      ###################################

      if (all(class(fit) != "try-error")){
        para<-as.data.frame(summary(fit)$coefficients)[-1,-3]
      }
      else{
        print(paste("Fitting problem for feature", x, "returning NA"))
        para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol = 3))
      }
      colnames(para)<-c('coef_LR', 'stderr_LR', 'pval_LR')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
    }
    return(para)
  })

  ###################
  # Combine results #
  ###################

  paras<-do.call(rbind, paras)

  #################
  # Return output #
  #################

  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)
}
