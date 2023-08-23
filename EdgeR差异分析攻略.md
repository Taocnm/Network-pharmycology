# EdgeR差异分析攻略

## step1:读取数据

## step2:构建分组矩阵，DGEList对象

```R
grouplist <- c(rep("sham",times = 3),rep("2VO",times = 3))
degs <- DGEList(counts = GSE199508_exprset, group = grouplist)
degs$samples
```

## step3：数据的预处理

1. 去除低表达量基因
2. 探索样本分组信息 -- 有助于挖掘潜在的差异样本

关于Normalization

在差异分析中，我们常常更关注的是相对表达量的变化，例如处理组的A基因表达量相对于对照组的而言是上调还是下调了。而基因表达量的定量准确性则在差异分析中不太重要，因此，在进行差异分析时，像RPKM/FPKM这种对转录本长度进行normalization方法是并不常用，也是没有必要的。

在常规的RNA-seq中，影响基因表达量更大的技术因素往往是**测序深度**以及**有效文库大小（effective libraries size）**。这也是一般的差异分析软件会进行normalize的部分。

```R
     # 1.过滤地表达
countsPerMillion <- cpm(degs)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
degs.keep <- degs[keep,]
dim(degs.keep)
    # 2.标准化
degs.norm <- calcNormFactors(degs.keep, method = 'TMM')
plotMDS(degs.norm, col=as.numeric(degs.norm$samples$group))
legend("bottomleft",as.character(unique(degs.norm$samples$group)), col=1:3, pch=20)
```

## step4:差异分析

```R
    # 1.构建分组矩阵
design <- model.matrix(~0+grouplist)
colnames(design) <- levels(factor(grouplist))
rownames(design) <- colnames(GSE199508_exprset)

    # 2.极大似然估计
degs.norm <- estimateGLMCommonDisp(degs.norm,design=design)
degs.norm <- estimateGLMTrendedDisp(degs.norm, design=design)
degs.norm <- estimateGLMTagwiseDisp(degs.norm, design=design)
plotBCV(degs.norm)

   # 3.统计检验
fit <- glmFit(degs.norm, design)
# LRT=likelihood ratio test
# group1-group2
lrt.1vs2 <- glmLRT(fit, contrast = c(1,-1))
degs.res.1vs3 <- topTags(lrt.1vs2, n = Inf, adjust.method = 'BH', sort.by = 'PValue')
DEG <- degs.res.1vs3$table

```

