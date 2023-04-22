    suppressMessages(library(stringr))
    suppressMessages(library(magrittr))
    suppressMessages(library(readr))
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))
    suppressMessages(library(dplyr))
    suppressMessages(library(Xmisc))
    suppressMessages(library(RNOmni))
    suppressMessages(library(cowplot))
    suppressMessages(library(psychometric))
    suppressMessages(library(bootstrap))
    source("/work-zfs/abattle4/ashton/prs_dev/prs_tools/liability_pseudoR2.R")

    getModels <- function(trait, n, data_, covars = "")
    {
        # function to obtain R-Squared from the data
        if (covars == "")
        {
            pval <- paste0('\`', n, '\`')
            models.lm <- c(as.formula(paste(trait, pval, sep = " ~ ")))
        } else { #we have covariates
            pname = paste0("\`",n,"\`")
            f_full <- as.formula(paste(trait, paste(pname,"+", covars), sep = " ~ "))
            f_part <-  as.formula(paste(trait, covars,sep = " ~ "))
            models.lm <- c(f_full, f_part)
        }
        return(models.lm)
    }
    iterFunct <- function(formulas, data, indices)
    {
            if (length(formulas) > 1)
            {
                lmv_complete <- lm(formulas[[1]], data =data[indices,])
                lmv_null <- lm(formulas[[2]], data =data[indices,])
                return(summary(lmv_complete)$r.squared - summary(lmv_null)$r.squared)
            }else{
                t <- lm(formulas[[1]], data =data[indices,])
                return(summary(t)$r.squared)
            }
     }
        


    jackknifeBootstrap <- function(trait, n, data, covars_)
    {
    models.lm <- getModels(trait, n, data, covars = covars_)
    j <- jackknife(1:nrow(data), iterFunct, data= data, formulas = models.lm)
    max_i <- which(j$jack.values == max(j$jack.values))
    min_i <- which(j$jack.values == min(j$jack.values))
    print(paste("Max R2: index = ", max_i, "value = ", j$jack.values[max_i]))
    print(data[max_i,1])
    print(paste("Diff from mean",j$jack.values[max_i] - mean(j$jack.values))) 
    #print(paste("Min R2: index = ", min_i, "value = ", j$jack.values[min_i]))
    #print(data[min_i,1])
    CI <- qt(1-0.025, nrow(data)-1) * j$jack.se
    return(c(CI, mean(j$jack.values)))
    }

    jackknife95 <- function(trait, n, data, covars_)
    {
    models.lm <- getModels(trait, n, data, covars = covars_)
    j <- jackknife(1:nrow(data), iterFunct, data= data, formulas = models.lm)
    max_i <- which(j$jack.values == max(j$jack.values))
    min_i <- which(j$jack.values == min(j$jack.values))
    CIs <- quantile(j$jack.values, probs = c(0.025,0.975))
    return(CIs)
    }


    #covars is the raw covars string
    straightBootstrap <- function(trait, n, data_, boot_type, covars)
    {
        #code based on https://www.statmethods.net/advstats/bootstrapping.html
        # Bootstrap 95% CI for R-Squared
        suppressMessages(library(boot))
        models.lm <- getModels(trait, n, data_, covars = covars)
        reps = nrow(data_) + 1
        #reps = 2000
        results <- boot(data = data_, statistic = iterFunct, R =reps, formulas = models.lm)
        max_i <- which(results$t == max(results$t))
        min_i <- which(results$t == min(results$t))
        print(paste("Max: index = ", max_i, "value = ", results$t[max_i]))
        print(data_[max_i,])
        print(paste("Min: index = ", min_i, "value = ", results$t[min_i]))
        print(data_[min_i,])

        # get 95% confidence interval
        t <-  boot.ci(results, type=boot_type)
        if (boot_type == "norm") {
            return(c(t$normal[2], t$normal[3]))
            } else if (boot_type == "bca"){
                return(c(t$bca[4], t$bca[5]))
                
            } else if (boot_type == "perc"){
                return(c(t$percent[4], t$percent[5]))
            } else if (boot_type == "basic") {
                return(c(t$basic[4], t$basic[5]))
            } else {
                rd <- t[4]
                return(c(rd[4], rd[5]))
            }
}

    getSNPCounts <- function(log_file)
    {
        alldat = readLines(log_file)
        t <- seq(from = length(alldat), to = 1, by = -1)
        counts <- c()
        pvals <- c()
        switch = FALSE
        for (i in t)
        {
        if (grepl("Number of final SNPs used at", alldat[i], fixed=  T))
        {
            pvals <- c(pvals, str_match(alldat[i], "at  ([[:digit:][:punct:]e]+) ")[2])
            counts <- c(counts, str_match(alldat[i], ": ([[:digit:]]+) $")[2])
            switch = TRUE
        }else{
            if(switch) #we just came out of the region, just the first one
            {break}
            switch = FALSE
        }
        }
        return(data.frame("pval_names" = as.factor(pvals), "snp_counts" = as.numeric(counts)))
    }


    plotCorr <- function(dat, output, style_name, category_var, no_pvals, per_snp=FALSE)
    {
        library("ggsci")
        dat$pval_names <- as.numeric(as.character(dat$pval_names))
        dat <- arrange(dat, pval_names)
        dat$pval_names <- as.factor(dat$pval_names)
        dat$logp <- -log10(dat$pval_beta)
        if(category_var == ""){
            base <- ggplot(dat, aes(x = pval_names, y = r2, fill =logp )) + geom_bar(stat = "identity") + 
                labs(fill = "-log10(pval)") + scale_fill_gsea()
            if(no_pvals) {base}
            else {base <- base + geom_text(aes(label=as.character(round(pval_beta, digits =3))), position=position_dodge(width = 0.9), vjust = -0.2)}
            
        } else {
            base <- ggplot(dat, aes(x = pval_names, y = r2, fill =category )) + geom_bar(stat = "identity", position = "dodge") + 
                scale_fill_discrete(name = category_var)
            if(no_pvals) {base <- base}
            else {base <- base  + geom_text(aes(label = as.character(round(pval_beta, digits = 3))), position = position_dodge(width = 0.9), vjust = -0.3, size = 3)}
        }
        if(!per_snp)
        { 
            base <- base + labs(x="P-value threshold", y = expression(paste(R ^ 2)))
        } else 
        {
            base <- base + labs(x="P-value threshold", y = expression(paste("Per-snp R" ^ "2")))
        }

        p <- base + theme_minimal_grid(5) + geom_errorbar(aes(x=pval_names, ymin=r2_bars_lower, ymax= r2_bars_upper)) + theme(text = element_text(size=20))
        if(!is.na(output))
        {
          ggsave(p, filename = paste0(output, "bar_plot.",style_name, ".png"), height = 7, width = 8.5) 
        }
          return(p)

    }
    #Make a quantile plot
    plotQuantile <- function(dat, trait, output, n_quants, style_name, covars_arg)
    {
        covars <- str_split(covars_arg, fixed("+"))[[1]]
        for (n in pval_names) #is pval_names global here?
        {
            for_plot <-  dat %>% mutate(quantile = ntile(dat[[n]], n_quants)) 
            if(covars[1] != "")
            {for_plot <- for_plot %>% dplyr::select(IID, all_of(n), quantile, all_of(trait), all_of(covars)) %>% filter(!(is.na(trait))) %>% ungroup()}
            else
            {for_plot <- for_plot %>% dplyr::select(IID, all_of(n), quantile, all_of(trait)) %>% filter(!(is.na(trait))) %>% ungroup()}
            xounts <- count(for_plot, quantile)
            avgs <- for_plot %>% group_by(quantile) %>% dplyr::summarize(avg = mean(!!as.symbol(trait)), std_dev = sd(!!as.symbol(trait))) %>% arrange(quantile) %>% mutate(num = xounts$n) %>% mutate(sem = std_dev/sqrt(num))
            #Simple one
            if (style_name == "simple")
            {
                plt <- (ggplot(avgs, aes(x=quantile, y=avg)) + geom_pointrange(aes(ymin=avg-sem, ymax=avg+sem)) + 
                            ylab(trait) + xlab("Quantile") + scale_x_continuous(breaks=c(1:n_quants), labels=c(1:n_quants)) + 
                            ggtitle(paste("Simple Quantile plot, p = ", n)) + theme_minimal_grid(15)) #+ scale_color_npg())
                
                ggsave(filename = paste0(output, ".quantile_plot.",n, ".", style_name, ".png"), plot = plt, height = 7, width = 8.5) 
            } else
            {
                #Nuanced one that allows for correcting for covariate- identify reference quantile- the median one.
                #med_val <- which(for_plot[[n]] == median(for_plot[[n]]))
                #med_quant <- for_plot[med_val,]$quantile
                med_quant <- n_quants/2#the one in the middle, right?
                beta <- c()
                std_error <- c()
                #$Compare each other quantile to the reference quantile
                for(i in 1:n_quants)
                {
                    dat_c <- for_plot %>% filter(quantile %in% c(i, med_quant)) %>% mutate(explan  = ifelse(quantile == i, 1,0))
                    #Still need to add lines for residual standard error, or stadnard error on beta, depending on what we choose.
                    if(i == med_quant)
                    {
                        beta <- c(beta, 0)
                        std_error <- c(std_error, 0)
                        next
                    }
                    if(covars[1] == "")
                    {
                        lmv_prs <- lm(dat_c[[trait]] ~ dat_c$explan)
                    }
                    else{
                        
                        f_covars <- as.formula(paste(trait, paste(covars), sep = " ~ "))
                        lmv_covars <- lm(f_covars, data = dat_c)
                        lmv_prs <- lm(lmv_covars$residuals ~ dat_c$explan) #Could also do it all together and just pull out the betas. Would be worth looking empirically if they are the same, which sthey should be GIVEN independence
                    }
                    beta <- c(beta, summary(lmv_prs)$coefficients[2,1])
                    std_error <- c(std_error, summary(lmv_prs)$coefficients[2,2])
                }
                quantile_covars <- data.frame(beta, std_error, quant = 1:n_quants) 
                plt <- (ggplot(quantile_covars, aes(x=quant, y=beta)) + geom_pointrange(aes(ymin=beta-(qnorm(1-0.025)*std_error), ymax=beta+(qnorm(1-0.025)*std_error))) + 
                            ylab("Effect size") + xlab("Quantile") +  scale_x_continuous(breaks=c(1:n_quants), labels=c(1:n_quants)) + 
                            ggtitle(paste("Relative Quantile plot, p = ", n, "Covariates:", covars_arg)) + theme_minimal_grid(15))
                ggsave(filename = paste0(output, ".quantile_plot.",n, ".", style_name, ".png"), plot = plt, height = 7, width = 8.5) 
            }
        
        }
    }

