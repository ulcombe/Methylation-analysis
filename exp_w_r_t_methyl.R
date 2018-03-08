# R module used to contain methods used for 
# testing expression(dasl) with respect to
# methylation status
# Author: love
###############################################################################
require(MASS)
require(prada)
require(gridExtra)
#
# Calculate mean methylation of tumour > 30% different from normal mean
# 	- if differential in only one direction classify
#
cpg_tumour_mean_diff <- function(x,phenoData,threshold){
	#x <- bval["cg00009292",] "cg17465631",] #cg27539986",] # x <- bval[1,]
	#threshold <- 0.3
	phenoData <- pd
	cpg_mean_tn <- tapply(x,
			as.character(phenoData$Sample_Group),
			mean)
	t_cpg <- x[which(as.character(phenoData$Sample_Group) == 'T')]
	diff_cpg <- t_cpg - cpg_mean_tn['N']
	diff_cpg_mean <- NA 
	cpg_calc_cat <- NA
	cpg_metcat <- NA
	cpg_t_95th <- quantile(x[which(phenoData$Sample_Group == 'T')], .95)
	cpg_t_5th <- quantile(x[which(phenoData$Sample_Group == 'T')], .05)
	
	if(length(which(diff_cpg > threshold)) > 0 & length(which(diff_cpg < -threshold)) > 0){
		diff_count <- NA
		
	}else if(length(which(abs(diff_cpg) > threshold)) ==0){
		diff_count <- length(which(abs(diff_cpg) > 0.3))		
	}else{
		diff_cpg_mean <- mean(diff_cpg[which(abs(diff_cpg) > 0.3)])
		cpg_mean_tn['T'] <- mean(t_cpg[which(abs(diff_cpg) > 0.3)])
		diff_count <- length(which(abs(diff_cpg) > 0.3))
		
		# classify
		cpg_calc_cat <- (cpg_mean_tn['T']) - (cpg_mean_tn['N'])
		cpg_metcat <- 'hypo'
		if(cpg_calc_cat > 0){
			cpg_metcat <- 'hyper'
		}
	}
	cat_df <- data.frame(
				n_mean=cpg_mean_tn[1],
				t_mean=cpg_mean_tn[2],
				t_95th=cpg_t_95th,
				t_5th=cpg_t_5th,
				calc_cat=cpg_calc_cat,
				category=as.character(cpg_metcat),
				mean_met_delta=diff_cpg_mean,
				count=diff_count)
	return(cat_df)
}

#
# calculate the methylation difference between mean normal value
# and mean tumour value of those samples >30% different from the 
# mean. Where the change methylation is only in one direction
#
cpg_mean_diff <- function(cpg,phenoData,threshold){
	#cpg <- "cg23244095" "cg15381475" "cg02279071" "cg17465631",] #cg27539986",] # x <- bval[1,]
	#threshold <- 0.3
	#phenoData <- pd
	x <- bval[cpg,]
	cpg_mean_tn <- tapply(x,
			as.character(phenoData$Sample_Group),
			mean)
	t_cpg <- x[which(as.character(phenoData$Sample_Group) == 'T')]
	diff_cpg <- t_cpg - cpg_mean_tn['N']
	diff_cpg_mean <- NA 
	cpg_calc_cat <- NA
	cpg_metcat <- NA
	cpg_t_95th <- quantile(x[which(phenoData$Sample_Group == 'T')], .95)
	cpg_t_5th <- quantile(x[which(phenoData$Sample_Group == 'T')], .05)
	normal_min <- min(x[which(phenoData$Sample_Group == 'N')])
	normal_max <- max(x[which(phenoData$Sample_Group == 'N')])
	cpg_mean_t_30pc <- NA
	exp_mean_t_30pc <- NA
	exp_mean_n <- NA
	expr_ttest_pval <- NA
	
	# Don't consider CpGs with beta > threshold(30%) in both hyper and hypo direction
	if(length(which(diff_cpg > threshold)) > 0 & length(which(diff_cpg < -threshold)) > 0){
		count_gt_threshold <- NA
		count_lt_threshold <- NA
	# Count number of samples greater than threshold  	
	}else if(length(which(abs(diff_cpg) > threshold)) == 0){
		count_gt_threshold <- length(which(abs(diff_cpg) > threshold))
		count_lt_threshold <- length(which(abs(diff_cpg) > threshold))
	}else{
		diff_cpg_mean <- mean(diff_cpg[which(abs(diff_cpg) > threshold)])
		cpg_mean_t_30pc <- mean(t_cpg[which(abs(diff_cpg) > threshold)])
		count_gt_threshold <- length(which(abs(diff_cpg) > threshold))
		count_lt_threshold <- length(which(abs(diff_cpg) < threshold))
		
		# Get expression data
		gene_id <- ruv_refseq_annot[match(cpg,ruv_refseq_annot[,'cpg.Name']),'cpg_island.txName']
		if(gene_id %in% rownames(dge_obj.rpkm.fil_T)){
			threshold_samples <- names(diff_cpg)[which(abs(diff_cpg) > threshold)]
			exp_t <- dge_obj.rpkm.fil_T[gene_id,which(colnames(dge_obj.rpkm.fil_T) %in% threshold_samples)]
			exp_n <- dge_obj.rpkm.fil_N[gene_id,]
			exp_mean_t_30pc <- mean(as.numeric(exp_t)) 
			exp_mean_n <- mean(as.numeric(exp_n))
			if(length(exp_t)> 1){
				expr_ttest_pval <- t.test(log2(exp_n+1),log2(exp_t+1))$p.value
			}else{
				expr_ttest_pval <- NA
			}
			
		}
		# classify
		cpg_calc_cat <- (cpg_mean_tn['T']) - (cpg_mean_tn['N'])
		cpg_metcat <- 'hypo'
		if(cpg_calc_cat > 0){
			cpg_metcat <- 'hyper'
		}
	}

	cat_df <- data.frame(
			n_mean=cpg_mean_tn[1],
			t_mean=cpg_mean_t_30pc,
			overall_t_mean=cpg_mean_tn['T'],
			t_95th=cpg_t_95th,
			t_5th=cpg_t_5th,
			normal_min,
			normal_max,
			exp_mean_t_30pc,
			exp_mean_n,
			expr_ttest_pval,
			calc_cat=cpg_calc_cat,
			category=as.character(cpg_metcat),
			mean_met_delta=diff_cpg_mean,
			count_lt_threshold,
			count_gt_threshold)

	return(cat_df)
}

categorise_mean <- function(cpg){
	# cpg <- "cg02279071" 
	# pd <- pData(pd)
	
	calc_bval_df <- cpg_mean_diff(cpg,pd,0.3)
	colnames(calc_bval_df) <- paste("b",colnames(calc_bval_df),sep="_")
	cat_df <- cbind(cpg=rownames(bval[cpg,,drop=F]),calc_bval_df)
	return(cat_df)
}


mval_foldchange <- function(x){
	
	#x <- 1 #'cg07073544;PRKG2'
	#x <- which(cat_test_all_df[,'combined_id'] == 'cg05048523;OBSCN') #'cg01302270;MLH1')
	#print(x['cpg'])
	
	class <-NA
	class_trim <- NA
	mval_fc_diff <- NA
	exp_fc_diff <- NA
	x <- cat_test_all_df[x,]
	if(x['b_category'] == 'hyper'){
		mval_fc_diff <- as.numeric(x['mval_fc_95th'])
		
		if(as.numeric(x['lm_coefficient']) > 0 & abs(as.numeric(x['exp_fc_95th'])) > abs(as.numeric(x['exp_fc_5th']))) {
			class <- 'Upreg_hyper_pos_pos'
			class_trim <- 'Upreg_hyper'
			
			exp_fc_diff <- as.numeric(x['exp_fc_95th'])
			
		}else if(as.numeric(x['lm_coefficient']) < 0 & abs(as.numeric(x['exp_fc_95th'])) < abs(as.numeric(x['exp_fc_5th']))){
			class <- 'Downreg_hyper_neg_pos'
			class_trim <- 'Downreg_hyper'
			exp_fc_diff <- (as.numeric(x['exp_fc_5th']))
			
		}else{
		
		}
		
	}else if(x['b_category'] == 'hypo'){
		
		mval_fc_diff <- (as.numeric(x['mval_fc_5th']))
		
		if(as.numeric(x['lm_coefficient']) > 0 & abs(as.numeric(x['exp_fc_95th'])) < abs(as.numeric(x['exp_fc_5th'])) ){
			class <- 'Downreg_hypo_neg_neg'
			class_trim <- 'Downreg_hypo'
			exp_fc_diff <- (as.numeric(x['exp_fc_5th']))
			
		}else if(as.numeric(x['lm_coefficient']) < 0 & abs(as.numeric(x['exp_fc_95th'])) > abs(as.numeric(x['exp_fc_5th'])) ){
			class <- 'Upreg_hypo_pos_neg'
			class_trim <- 'Upreg_hypo'
			exp_fc_diff <- as.numeric(x['exp_fc_95th'])
			
		}else{
		
		}
		
	}else{
		print("Error")
	}
	mode(class_trim) <- 'character'
	mode(class) <- 'character'
	fc_df <- data.frame(mval_fc_diff,exp_fc_diff,class,class_trim)
	return(fc_df) 
}

#
# count number in each quadrant 
#
count_quadrant <- function(x,met_col,dasl_col){
	count_table <-  thresholds(as.numeric(x[,met_col]),
			as.numeric(x[,dasl_col]),xthr=0,ythr=0)
	prop_table <- round(prop.table(count_table)*100,digits=1)
	pval <- round(fisher.test(count_table)$p.value,digits=3)
	rownames(count_table) <- c('exp_up','exp_down')
	rownames(prop_table) <- c('prop_exp_up','prop_exp_down')
	colnames(count_table) <- c('met_hypo','met_hyper')
	colnames(prop_table) <- c('met_hypo','met_hyper')
	return(list(count_table,prop_table,pval))
}

#
# linear models for tcga validation
#
lm_mval_met_exp_tcga_validate <- function(row_x,figure){
	
	#row_x <- cpg2gene_intersect_df[1,]"cg14751544;MLH1",,drop=T]
	cpg_probe <- as.character(row_x[['cpg']])
	exp_symbol <- as.character(row_x[['cpg_island.txName']])
	#cpg_probe <- 'cg09296001'
	#exp_symbol <- 'SND1'
	
	# get met change v matched normal samples
	int_df <- data.frame(
			mval=as.numeric(tumour_mval[cpg_probe,,drop=T]),
			bval=as.numeric(tumour_bval[cpg_probe,,drop=T]),
			exp_symbol=as.numeric(exp_log2_norm_matched[exp_symbol,,drop=T]),
			exp_logfc=as.numeric(diff_exp_df[exp_symbol,,drop=T]),
			mval_fc=as.numeric(t_m_diff_matched_df[cpg_probe,,drop=T]),
			class='NA')
	
	colnames(int_df) <- c('mval','bval','exp_symbol','exp_logfc','mval_fc','class')
	rownames(int_df) <- colnames(tumour_mval)
	int_df$class <- factor('N-like',levels=c('N-like','T-like'))
	
	#linear model of all samples dasl v met
	fit <- lm(int_df[,'exp_logfc'] ~ int_df[,'mval_fc'] )
	rfit <- rlm(int_df[,'exp_logfc'] ~ int_df[,'mval_fc'] )
	conf_25 <- confint(fit)[1,1]
	conf_975 <- confint(fit)[1,2]
	
	#test mean of tumours based on methylation status
	if(row_x$b_category == 'hyper'){
		exp_t_change <- int_df[which(int_df$bval > row_x$b_normal_max),'exp_symbol']
		exp_t_norm <- int_df[which(int_df$bval < row_x$b_normal_max),'exp_symbol']
		int_df$class[which(int_df$bval > row_x$b_normal_max)] <- 'T-like'
	}else{
		exp_t_change <- int_df[which(int_df$bval < row_x$b_normal_min),'exp_symbol']
		exp_t_norm <- int_df[which(int_df$bval > row_x$b_normal_min),'exp_symbol']
		int_df$class[which(int_df$bval < row_x$b_normal_max)] <- 'T-like'
	}
	ttest_pval <- NA
	
	if(length(exp_t_change) > 1 & length(exp_t_norm) > 1){
		ttest_pval <- t.test(exp_t_change,exp_t_norm)$p.value
	}
	#spearmans rank
	sp_test <- cor.test(int_df$mval,int_df$exp_symbol,method = "spearman")
	
	# get 5th and 95th quantile ranges
	mval_q <- quantile(int_df$mval_fc, c(.05, .95))
	exp_q <- quantile(int_df$exp_logfc, c(.05, .95))
	
	if(figure){
		
		p1 <- ggplot(int_df, aes(x=mval_fc, y=exp_logfc,col=class)) +
				geom_point(size=3)  +
				geom_abline(intercept = fit$coefficient[1], slope = fit$coefficient[2],col="red") +
				geom_abline(intercept = rfit$coefficient[1], slope = rfit$coefficient[2],col="red",linetype="dotted") +
				#geom_abline(intercept = rfit$coefficient[1], slope = rfit$coefficient[2],col="goldenrod") +
				#geom_point(data=data.frame(beta=beta_q,dasl=dasl_q), aes(x=beta, y=dasl), colour="red", size=5) +
				#scale_y_continuous(limits=c(-6, 6)) +
				#scale_x_continuous(limits=c(-7, 7)) +
				#geom_density2d() +
				geom_hline(yintercept=0) +
				geom_vline(xintercept=0) +
				theme_bw()
		p2 <- ggplot(int_df, aes(x=bval, y=exp_symbol,col=class)) +
				geom_point(size=3)  +
				#geom_abline(intercept = rfit$coefficient[1], slope = rfit$coefficient[2],col="goldenrod") +
				#geom_point(data=data.frame(beta=beta_q,dasl=dasl_q), aes(x=beta, y=dasl), colour="red", size=5) +
				#scale_y_continuous(limits=c(-6, 6)) +
				#scale_x_continuous(limits=c(-7, 7)) +
				#geom_density2d() +
				geom_hline(yintercept=0) +
				geom_vline(xintercept=0) +
				theme_bw()
		print(grid.arrange(p1,p2,ncol=2))
		#png(paste('tcga_',cpg_probe,'_',exp_symbol,'_lm.png',sep=""))
		#print(p1)
		#dev.off()
		#abline(fit,col="red")
		#dev.off()
		#X11();
		#layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
		#plot(fit)							
	}
	#summary(fit)
	
	# meth_mean_norm_nc <- mean(nc_met_row_pair)
	# meth_mean_norm_c <- mean(c_met_row_pair)
	
	rlm_p_value =  2*pt(abs(summary(rfit)$coefficients[2,'t value']), summary(rfit)$df[2], lower.tail=FALSE) 
	
	lm_df <- data.frame(cpg=cpg_probe,
			exp_symbol=exp_symbol,
			r_squared=summary(fit)$r.squared,
			lm_coefficient=unlist(fit$coefficients[2]),
			lm_p_value=anova(fit)$"Pr(>F)"[1],
			lm_conf_2_5pc=conf_25,
			lm_conf_97_5pc=conf_975,
			rlm_p_value, 
			rho=sp_test$estimate,
			sp_p_value=sp_test$p.value,
			mval_fc_5th=mval_q[1],
			mval_fc_95th=mval_q[2],
			exp_fc_5th=exp_q[1],
			exp_fc_95th=exp_q[2],
			ttest_pval,
			mean_exp_t_change=mean(exp_t_change),
			mean_exp_t_norm=mean(exp_t_norm)	
	)
	return(lm_df)
}


#
# linear models for cell line validation
#
lm_mval_met_exp_5AZA <- function(cpg_probe,exp_symbol,figure){
	
	#cpg_probe <- 'cg05925365'
	#exp_symbol <- 'ICMT'
	
	# get met change v matched normal samples
	int_df <- data.frame(logfc_5aza_exp=as.double(logfc_5AZA_df[exp_symbol,,drop=T]),
						diff_mval=as.double(diff_MSet[cpg_probe,,drop=T]))
	
	#linear model of all samples dasl v met
	fit <- lm(data=int_df, logfc_5aza_exp ~ diff_mval)
	conf_25 <- confint(fit)[1,1]
	conf_975 <- confint(fit)[1,2]
	
	#spearmans rank
	sp_test <- cor.test(int_df$diff_mval,int_df$logfc_5aza_exp,method = "spearman")
	
	# get 5th and 95th quantile ranges
	mval_q <- quantile(int_df$diff_mval, c(.05, .95))
	exp_q <- quantile(int_df$logfc_5aza, c(.05, .95))
	
	if(figure){
		
		p1 <- ggplot(int_df, aes(x=diff_mval, y=logfc_5aza_exp)) +
				geom_point(size=3,shape=1)  +
				geom_abline(intercept = fit$coefficient[1], slope = fit$coefficient[2],col="red") +
				#geom_point(data=data.frame(beta=beta_q,dasl=dasl_q), aes(x=beta, y=dasl), colour="red", size=5) +
				#scale_y_continuous(limits=c(-6, 6)) +
				#scale_x_continuous(limits=c(-7, 7)) +
				geom_hline(yintercept=0) +
				geom_vline(xintercept=0) +
				theme_bw()
		print(p1)
		png(paste(config$concordant$outdir,'/5aza_',cpg_probe,'_',exp_symbol,'_lm.png',sep=""))
		print(p1)
		dev.off()
		#abline(fit,col="red")
		#dev.off()
		#X11();
		#layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
		#plot(fit)							
	}
	summary(fit)
	
	# meth_mean_norm_nc <- mean(nc_met_row_pair)
	# meth_mean_norm_c <- mean(c_met_row_pair)
	
	lm_df <- data.frame(cpg=cpg_probe,
			exp_symbol=exp_symbol,
			r_squared=summary(fit)$r.squared,
			lm_coefficient=unlist(fit$coefficients[2]),
			lm_p_value=anova(fit)$"Pr(>F)"[1],
			lm_conf_2_5pc=conf_25,
			lm_conf_97_5pc=conf_975,
			rho=sp_test$estimate,
			sp_p_value=sp_test$p.value,
			mval_5th=mval_q[1],
			mval_95th=mval_q[2],
			exp_logfc_5th=exp_q[1],
			exp_logfc_95th=exp_q[2]
	)
	return(lm_df)
}
