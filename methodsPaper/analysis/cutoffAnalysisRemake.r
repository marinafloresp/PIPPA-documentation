library(kneedle)
library(segmented)
library(tidyr)
library(dplyr)
options(scipen = 5)

peakDir <-"/dir/to/peak/set/"
peakFile <- paste0(peakDir, "/", "peakSetName")

## Load peak set and add column names
peaks.gr<-makeGRangesFromDataFrame(read.table(file.path(peakFile)), keep.extra.columns=TRUE, seqnames.field="V1",
                                   start.field="V2",
                                   end.field="V3",
                                   ignore.strand=TRUE)
peaks.df <- as.data.frame(peaks.gr)
colnames(peaks.df) <- c("chr", "start", "end", "width", "strand", "peakID", "score", "strand2", "signal", "pvalue", "qvalue", "score2")
peaks.df$pvals <- 10^(-peaks.df$pvalue)
peaks.df$qvals<- 10^(-peaks.df$qvalue)

## Set of thresholds to be considered. From 0.05 to 1e-100
exponents <- seq(-2, -100, length.out = 50)
thresholds <- c(0.05,10^exponents)
order_levels.og <- as.character(thresholds)

## Number of peaks vs q-values
qvals <- peaks.df$qvals
counts_list <-  sapply(thresholds, function(t) sum(qvals <= t, na.rm = TRUE))

df <- data.frame(thres = thresholds,count = counts_list)
numPeaks <- df$count
qvals <- -log10(df$thres)

## calculate elbow
knee <-kneedle(qvals, df$count)
knee_log <-  format(10^-knee[1], scientific = TRUE)

## segmented analysis: fits two linear regimes and estimates a breakpoint.
model <- lm(numPeaks ~ qvals)
seg_model <- segmented(model, seg.Z = ~qvals)
summary(seg_model)
seg_out<- round(abs(seg_model$psi[2]))

ci_segmented <- confint(seg_model)
ci_segmented1 <- format(10^-round(abs(ci_segmented[1])), scientific = TRUE)
ci_segmented2 <- format(10^-round(abs(ci_segmented[2])), scientific = TRUE)
ci_segmented3 <- format(10^-round(abs(ci_segmented[3])), scientific = TRUE)

## Plot 1: Number of peaks vs q-values
order_levels <- c(order_levels.og,ci_segmented1,ci_segmented2,ci_segmented3)
order_levels <- as.factor( sort(as.numeric(order_levels), decreasing = TRUE))

p1<- ggplot(df, aes(x = as.factor(thres), y = count, group=1)) +
  geom_point() +
  geom_line(size = 1)+
  labs(
    x = "q-value",
    y = "Number of peaks", title = "Cutoff analysis (Number of peaks)"
  )+
  geom_vline(color='blue',size = 1.5,aes(xintercept = as.factor(knee_log)))+
  geom_vline(color='red',size = 1.5,aes(xintercept = as.factor(ci_segmented1)))+
  geom_vline(color='#FF6A6A',size = 1.5,aes(xintercept = as.factor(ci_segmented2)))+
  geom_vline(color='#FF6A6A',size = 1.5,aes(xintercept = as.factor(ci_segmented3)))+
  theme_minimal(base_size = 14) +scale_x_discrete(limits = order_levels)+
  theme(
    axis.text.x = element_text(size = 10,angle = 90,face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold")
  )

## Total length of peaks vs q-values
df.plot2 <- tibble(thres = thresholds) %>%
                       rowwise() %>%
  mutate(total_width = sum(peaks.df$width[peaks.df$qvals <= threshold])) %>% 
                       as.data.frame()
total_width <- df.plot2$total_width
qvals <- -log10(df.plot2$thres)

## calculate elbow
knee <-kneedle(qvals, df.plot2$total_width)
knee_log <- format(10^-knee[1], scientific = TRUE)

## segmented analysis
model <- lm(total_width ~ qvals)
seg_model <- segmented(model, seg.Z = ~qvals)
summary(seg_model)
seg_out<- round(abs(seg_model$psi[2]))
ci_segmented <- confint(seg_model)
ci_segmented1 <- format(10^-round(abs(ci_segmented[1])), scientific = TRUE)
ci_segmented2 <- format(10^-round(abs(ci_segmented[2])), scientific = TRUE)
ci_segmented3 <- format(10^-round(abs(ci_segmented[3])), scientific = TRUE)

## Plot 2: Total length of peaks vs q-values
order_levels <- c(order_levels.og,ci_segmented1,ci_segmented2,ci_segmented3)
order_levels <- as.factor( sort(as.numeric(order_levels), decreasing = TRUE))

p2<-ggplot(df.plot2, aes(x = factor(threshold), y = total_width, group=1)) +
  geom_point() +
  geom_line(size = 1)+
  labs(
    x = "q-value",
    y = "Total length peaks (bp)", title ="Cutoff analysis (Total length peaks)"
  )+
  geom_vline(color='blue',size = 1.5,aes(xintercept = as.factor(knee_log)))+
  geom_vline(color='red',size = 1.5,aes(xintercept = as.factor(ci_segmented1)))+
  geom_vline(color='#FF6A6A',size = 1.5,aes(xintercept = as.factor(ci_segmented2)))+
  geom_vline(color='#FF6A6A',size = 1.5,aes(xintercept = as.factor(ci_segmented3)))+
  theme_minimal(base_size = 14) +scale_x_discrete(limits = order_levels)+
  theme(
    axis.text.x = element_text(size = 10,angle = 90,face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold")
  )
                       
## Average length of peaks vs q-values
df.plot3 <- tibble(thres = thresholds) %>%
  rowwise() %>%
  mutate(mean_width = round(mean(peaks.df$width[peaks.df$qvals <= thres])))

mean_width <- df.plot3$mean_width
qvals <- -log10(df.plot3$thres)
 
## calculate elbow
knee <-kneedle(qvals, mean_width, decreasing = FALSE)
knee_log <- format(10^-knee[1], scientific = TRUE)

## segmented analysis
model <- lm(mean_width ~ qvals)
seg_model <- segmented(model, seg.Z = ~qvals)
summary(seg_model)
seg_out<- round(abs(seg_model$psi[2]))
ci_segmented <- confint(seg_model)
ci_segmented1 <- format(10^-round(abs(ci_segmented[1])), scientific = TRUE)
ci_segmented2 <- format(10^-round(abs(ci_segmented[2])), scientific = TRUE)
ci_segmented3 <- format(10^-round(abs(ci_segmented[3])), scientific = TRUE)

## Plot 3: Average length of peaks vs q-values
order_levels <- c(order_levels.og,ci_segmented1,ci_segmented2,ci_segmented3)
order_levels <- as.factor( sort(as.numeric(order_levels), decreasing = TRUE))

p3 <-ggplot(df.plot3, aes(x = factor(thres), y = mean_width, group=1)) +
  geom_point() +
  geom_line(size = 1)+
  labs(
    x = "q-value",
    y = "Average length peaks (bp)", title ="Cutoff analysis (mean length peaks)"
  )+
  geom_vline(color='blue',size = 1.5,aes(xintercept = as.factor(knee_log)))+
  geom_vline(color='red',size = 1.5,aes(xintercept = as.factor(ci_segmented1)))+
  geom_vline(color='#FF6A6A',size = 1.5,aes(xintercept = as.factor(ci_segmented2)))+
  geom_vline(color='#FF6A6A',size = 1.5,aes(xintercept = as.factor(ci_segmented3)))+
  theme_minimal(base_size = 14) +scale_x_discrete(limits = order_levels)+
  theme(
    axis.text.x = element_text(size = 10,angle = 90,face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold")
  )
