library(DEP)
set.seed(12345)
MY_test_diff_in_DEP <- function (se, type = c("control", "all", "manual"), control = NULL, 
                                 test = NULL, design_formula = formula(~0 + condition)) 
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), 
                          is.character(type), class(design_formula) == "formula")
  type <- match.arg(type)
  col_data <- colData(se)
  raw <- assay(se)
  if (!is.null(control)) {
    assertthat::assert_that(is.character(control), length(control) == 
                              1)
  }
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>% 
    as.character() %>% .[-1]
  
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))

  if (type == "control") {
    if (is.null(control)) 
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control], 
                    control, sep = " - ")
  }
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_", 
                                           cntrst), collapse = ", "))
  
  fit            <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit   <- contrasts.fit(fit, made_contrasts)
  
  eB_fit <- eBayes(contrast_fit)
  
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf, 
                    confint = TRUE)
    res <- res[!is.na(res$t), ]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }
  cat(eB_fit$df.total)
  limma_res <- map_df(cntrst, retrieve_fun)
  table <- limma_res %>% select(rowname, logFC, CI.L, CI.R, 
                                P.Value, qval, comparison) %>% mutate(comparison = gsub(" - ", 
                               "_vs_", comparison)) %>% gather(variable, value, -c(rowname, 
                               comparison)) %>% mutate(variable = recode(variable, 
                               logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>% 
    unite(temp, comparison, variable) %>% spread(temp, value)
  rowData(se) <- merge(rowData(se), table, by.x = "name", 
                       by.y = "rowname", all.x = TRUE)
  return(list(se,eB_fit))
}

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

require(SummarizedExperiment);require(limma);require(purrr);require(tidyverse)
#BiocInstaller::biocLite("DEP")
library("DEP")
citation("DEP")

# For LFQ analysis
#run_app("LFQ")

# Loading a package required for data handling
library("dplyr")

data(UbiLength)
# The data is provided with the package
data <- UbiLength

# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively. 
data <- filter(data, Reverse != "+", Potential.contaminant != "+")
data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Are there any duplicated names?
data$name %>% duplicated() %>% any()

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Let's have a look at the SummarizedExperiment object
data_se


# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)



# Normalize the data
data_norm <- normalize_vsn(data_filt)



# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)


# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")

## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl, Ubi1_vs_Ctrl

as_tibble(data_diff@elementMetadata)



# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

## Tested contrasts: Ubi4_vs_Ubi6, Ubi4_vs_Ctrl, Ubi4_vs_Ubi1, Ubi6_vs_Ctrl, Ubi6_vs_Ubi1, Ctrl_vs_Ubi1

plot(data_diff_all_contrasts@elementMetadata[,c("Ctrl_vs_Ubi1_p.adj")])
sum(data_diff_all_contrasts@elementMetadata[,c("Ctrl_vs_Ubi1_p.adj")]<.05)


# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))

## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl

#Finally, significant proteins are defined by user-defined cutoffs using add_rejections.

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))





#################################
MYdata_diff <- MY_test_diff_in_DEP(se = data_imp, 
                                   type = "control", 
                                   control = "Ctrl")

XXX <- MYdata_diff[[2]]


df.total <- 10.29881

lZ <- qnorm(pt(XXX$t,df = XXX$df.total,log.p = T),log.p = T)

UBI1 <- lZ[,3]
UBI4 <- lZ[,1]
UBI6 <- lZ[,2]

hist(UBI1,breaks = 25,freq=F);curve(dnorm(x),add=T)
hist(UBI4,breaks = 25,freq=F);curve(dnorm(x),add=T)
hist(UBI6,breaks = 25,freq=F);curve(dnorm(x),add=T)

#################################################
#################################################
#################################################
pUBI1 <- 2*pnorm(abs(UBI1),lower.tail = F)
hist(pUBI1)
hist(p.adjust(pUBI1,"BH"))
sum(p.adjust(pUBI1,"BH")<.05)

pUBI4 <- 2*pnorm(abs(UBI4),lower.tail = F)
hist(pUBI4)
hist(p.adjust(pUBI4,"BH"))
sum(p.adjust(pUBI4,"BH")<.05)

pUBI6 <- 2*pnorm(abs(UBI6),lower.tail = F)
hist(pUBI6)
hist(p.adjust(pUBI6,"BH"))
sum(p.adjust(pUBI6,"BH")<.05)

#################################################
#################################################
#################################################
#################################################

Rcpp::sourceCpp('~/Git_repo_Ahead/Nollik2GM/Grouped_model/Hier_W3_v3.cpp')    
source('~/Git_repo_Ahead/Nollik2GM/Grouped_model/Hier_W3_v3.R')    

zz <- c(UBI1,UBI4,UBI6)
zg <- rep(1:3,rep(1899,3))

HP <- list(
  KAPPA = 2,
  aa    = 20,    bb   = 57,
  a_rho = 1,    b_rho = 9,
  ax    = 1,     bx   = 1, 
  m01   = -3,    m02  = 3,  # meno 5 5 è solution percorribile?
  V01   = 1,     V02  = 1, 
  a01   = 2,     b01  = 5,
  a02   = 2,     b02  = 5,
  a_phi0 = 10, b_phi0 = 10, 
  m_phi0 = 0,  V_phi0 = 1/100
)

G <- length(unique(zg))
res1 <- Nollik_w3_hier_v3(z  = zz, 
                          zg = zg,  
                          hyperpar = HP, 
                          NSIM = 10000,
                          BURN = 100000,
                          THIN = 50,
                          verbose.step = 50,
                          SIG1 = array(.5*diag(2),c(2,2,G)), 
                          SIG2 = array(.5*diag(2),c(2,2,G)),  
                          sigA = rep(.5,G))
###############################################################################################
saveRDS(object = res1,file = "~/Git_repo_Ahead/Nollik2GM/Grouped_model/New_limma_pm3_th50_BI100k.RDS")
#saveRDS(object = res1,file = "~/Git_repo_Ahead/Nollik2GM/Grouped_model/New_limma_pm5.RDS")

colMeans(res1$THcon)
plot(res1$RHOcon)
z <- zz
plot.prob <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  ef  <- locfdr::locfdr(z[zg==g],plot = F)
  pef <- 1-ef$fdr
  p1  <- apply(res1$Lcon[zg==g,],1,mean)
  points(p1~z[zg==g],col=2)  
  points(pef~z[zg==g],col=4)
}

plot.prob(1)
plot.prob(2)
plot.prob(3)

# plot(res1$M1con)
# plot(res1$M2con)
# plot(res1$M0)
# plot(res1$a)

##############################################################################################################
wd_mixmom_theta <- function(z,theta,m1,m2,k,s1=1,s2=1,a=1, logscale=F){
  g <- function(z)   {
    (1-exp(-(z/a)^(2*k)))*
      ( (1-theta) * dnorm(z,m1,sqrt(s1)) + theta * dnorm(z,m2,sqrt(s2)) )
  }
  Kost <- integrate(g,-40,40)$v
  RR <-   log(g(z)) - log(Kost)
  if(logscale){
    return(RR)
  }else{
    return(exp(RR))
  }
}
################################################################################################################
f_postj <- function(x,res1,j){
  (1-mean(res1$RHOcon)) * dnorm(x,mean(res1$M0[,j]),sqrt(mean(res1$S0[,j]))) +  
    mean(res1$RHOcon)  * wd_mixmom_theta(z = x,theta = mean(res1$THcon[,j]),
                                         m1 =    mean(res1$M1con[,j]),
                                         m2 =    mean(res1$M2con[,j]),
                                         s1 =    mean(res1$S1con[,j]),
                                         s2 =    mean(res1$S2con[,j]),
                                         k  =    round(mean(res1$Kcon)),
                                         a  =    mean(res1$a[,j]))
}
################################################################################################################

x <- seq(-7,7,length.out = 1000)
f1 <- f_postj(x,res1,1)
f2 <- f_postj(x,res1,2)
f3 <- f_postj(x,res1,3)
EFF <- list(f1,f2,f3)

plot.dens <- function(g){
  hist(z[zg==g],breaks = 40,freq=F,ylim=c(0,1))
  points(EFF[[g]]~x,col=2,pch=".")  
}
plot.dens(1)
plot.dens(2)
plot.dens(3)


f1_postj <- function(x,res1,j){
  mean(res1$RHOcon)  * wd_mixmom_theta(z = x,theta = mean(res1$THcon[,j]),
                                       m1 =    mean(res1$M1con[,j]),
                                       m2 =    mean(res1$M2con[,j]),
                                       s1 =    mean(res1$S1con[,j]),
                                       s2 =    mean(res1$S2con[,j]),
                                       k  =    round(mean(res1$Kcon)),
                                       a  =    mean(res1$a[,j]))
}
f0_postj <- function(x,res1,j){
  (1-mean(res1$RHOcon)) * dnorm(x,mean(res1$M0[,j]),sqrt(mean(res1$S0[,j]))) 
}

EFFE1 <- list()
EFFE0 <- list()

for(g in 1:G){
  
  EFFE0[[g]] <- f0_postj(x,res1,g)
  EFFE1[[g]] <- f1_postj(x,res1,g)
}



plot.dens2 <- function(g){
  hist(z[zg==g],breaks = 40,freq=F)
  points(EFFE0[[g]]~x,col=2,pch=".")  
  points(EFFE1[[g]]~x,col=4,pch=".")  
}
plot.dens2(1)
plot.dens2(2)
plot.dens2(3)

