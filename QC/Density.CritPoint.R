# Density estimation (automatically chosen bandwidth)
kdde_0 <- ks::kdde(x = pa4$nCount_RNA, deriv.order = 0)
kdde_0 <- dat5a.frame(estimate = kdde_0[["estimate"]],
                     eval.points = kdde_0[["eval.points"]])

# Density derivative estimation (automatically chosen bandwidth, but different
# from kdde_0!)
kdde_1 <- ks::kdde(x = pa4$nCount_RNA, deriv.order = 1)
kdde_1 <- data.frame(estimate = kdde_1[["estimate"]],
                     eval.points = kdde_1[["eval.points"]])

# Find point to place cut-off between empty droplets and cells
gradient_sign <- rle(kdde_1[["estimate"]]>0)
nf_cutoff <- kdde_1[["eval.points"]][sum(gradient_sign[["lengths"]][1:2])]

ed4.3 <- identify_empty_drops(nf_umi=test4, include_plot = TRUE, umi_rescue = 5700) # umi_rescue determined by local minima of nCount_RNA density distribution

