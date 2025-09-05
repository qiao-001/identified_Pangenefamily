# 加载所需包

library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(gridExtra)
setwd("d:/13.pan_genome/5.Camelina_sativa_LBDpan/LBD/14KaKS/")
# 读取数据文件
data <- read.csv("Kaksdata.txt",header = T,sep = "\t")
###eRNA和lncRNA的Cor的差异

p_Ka <-ggplot(data,aes(x=Coregene,y=Ka))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,
               aes(fill=Coregene),
               show.legend = F)+
  scale_fill_manual(values = c("#5B859E", "#1E395F", "#DF8D71","#75884B", "#1E5A46", "#DF8D71", "#AF4F2F", "#D46F90"))+
  theme_bw()+xlab("") + ylab(expression("Ka"))+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(axis.text = element_text(size = 14),axis.title=element_text(size=16))+
  theme(panel.grid=element_blank(),panel.border = element_rect(linewidth = 0.9))
p_Ka

p_Ka1 <-p_Ka+geom_signif(comparisons = list(c("Core", "SoftCore"), c("Core", "Dispensable"), c("SoftCore", "Dispensable")),
                                     y_position = c(0.09,0.10,0.08),
                                     vjust =0, #指定标记中文字部分与横线之间的距离
                                     test = wilcox.test,
                                     map_signif_level = function(p) sprintf("p = %.3e", p))
p_Ka1
ggsave("p_Ka1.tif",plot = p_Ka1,width = 800,height = 1000,units = "px")

p_Ks <-ggplot(data,aes(x=Coregene,y=Ks))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,
               aes(fill=Coregene),
               show.legend = F)+
  scale_fill_manual(values = c("#5B859E", "#1E395F", "#DF8D71","#75884B", "#1E5A46", "#DF8D71", "#AF4F2F", "#D46F90"))+
  theme_bw()+xlab("") + ylab(expression("Ks"))+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(axis.text = element_text(size = 14),axis.title=element_text(size=16))+
  theme(panel.grid=element_blank(),panel.border = element_rect(linewidth = 0.9))
p_Ks

p_Ks1 <-p_Ks+geom_signif(comparisons = list(c("Core", "SoftCore"), c("Core", "Dispensable"), c("SoftCore", "Dispensable")),
                         y_position = c(0.20,0.22,0.24),
                         vjust = 0, #指定标记中文字部分与横线之间的距离
                         test = wilcox.test,
                         map_signif_level = function(p) sprintf("p = %.3e", p))
p_Ks1
ggsave("p_Ks1.tif",plot = p_Ks1,width = 800,height = 1000,units = "px")




p_Ka_Ks <-ggplot(data,aes(x=Coregene,y=Ka_Ks))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(outlier.shape = 1,
               aes(fill=Coregene),
               show.legend = F)+
  scale_fill_manual(values = c("#5B859E", "#1E395F", "#DF8D71","#75884B", "#1E5A46", "#DF8D71", "#AF4F2F", "#D46F90"))+
  theme_bw()+xlab("") + ylab(expression("Ka/Ks"))+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(axis.text = element_text(size = 14),axis.title=element_text(size=16))+
  theme(panel.grid=element_blank(),panel.border = element_rect(linewidth = 0.9))

p_Ka_Ks

p_Ka_Ks1 <-p_Ka_Ks+geom_signif(comparisons = list(c("Core", "SoftCore"), c("Core", "Dispensable"), c("SoftCore", "Dispensable")),
                         y_position = c(2.7,3.05,3.35),
                         vjust = 0, #指定标记中文字部分与横线之间的距离
                         test = wilcox.test,
                         map_signif_level = function(p) sprintf("p = %.3e", p))
p_Ka_Ks1
ggsave("p_Ka_Ks1.tif",plot = p_Ka_Ks1,width = 800,height = 1000,units = "px")


allkaks<-grid.arrange(p_Ka1, p_Ks1, p_Ka_Ks1, ncol = 3)

ggsave("allkaks.tif",plot = allkaks,width = 2400,height = 1500,units = "px")





library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 模拟数据（替换为实际数据）

# 2. 计算Y轴范围和转换系数（适配双轴）
ka_range <- range(data$Ka)
ks_range <- range(data$Ks)
ka_buffer <- diff(ka_range) * 0.1
ks_buffer <- diff(ks_range) * 0.1

# 主坐标轴（Ks）范围
main_y_range <- c(ks_range[1] - ks_buffer, ks_range[2] + ks_buffer)
# 转换系数：使Ka适配主坐标轴范围
trans_coef <- max(main_y_range) / max(ka_range + ka_buffer)

# 3. 计算分组均值
stats_data <- data %>%
  group_by(Class) %>%
  summarise(
    Ka_mean = mean(Ka, na.rm = TRUE),
    Ks_mean = mean(Ks, na.rm = TRUE),
    KaKs_mean = mean(Ka_Ks, na.rm = TRUE),
    KaKs_sd = sd(Ka_Ks, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    KaKs_lower = KaKs_mean - KaKs_sd,
    KaKs_upper = KaKs_mean + KaKs_sd,
    Ka_mean_scaled = Ka_mean * trans_coef  # 缩放Ka均值以适配主坐标轴
  )

# 4. 处理Ka/Ks数据（调整高度适配主轴）
kaks_data <- stats_data %>%
  select(Class, KaKs_mean, KaKs_lower, KaKs_upper) %>%
  mutate(
    KaKs_mean_adjust = KaKs_mean * (max(main_y_range) / max(KaKs_mean * 1.2))
  )

# 5. 确保Class顺序
all_classes <- unique(data$Class)
stats_data$Class <- factor(stats_data$Class, levels = all_classes)
kaks_data$Class <- factor(kaks_data$Class, levels = all_classes)

# 6. 绘制图形（Ka和Ks左右分布，不加点）
# 续接前文数据处理部分...

# 绘制图形
linesaa<-ggplot() +
  # 柱状图（左右分布）
  geom_col(
    data = stats_data,
    aes(x = Class, y = Ka_mean_scaled, fill = "Ka"),
    width = 0.35,
    position = position_nudge(x = -0.2),
    linewidth = 0.3,
    color = "black",
    alpha = 0.85
  ) +
  geom_col(
    data = stats_data,
    aes(x = Class, y = Ks_mean, fill = "Ks"),
    width = 0.35,
    position = position_nudge(x = 0.2),
    linewidth = 0.3,
    color = "black",
    alpha = 0.85
  ) +
  
  # Ka/Ks趋势区域（带fill映射，用于图例）
  geom_ribbon(
    data = kaks_data,
    aes(x = Class, ymin = KaKs_lower, ymax = KaKs_upper, fill = "Ka/Ks"),
    alpha = 0.25,
    color = NA
  ) +
  
  # Ka/Ks折线（带color映射，用于图例）
  geom_line(
    data = kaks_data,
    aes(x = Class, y = KaKs_mean_adjust, group = 1, color = "Ka/Ks"),
    size = 1.0,
    lineend = "round"
  ) +
  
  # Ka/Ks数据点和标签
  geom_point(
    data = kaks_data,
    aes(x = Class, y = KaKs_mean_adjust),
    color = "red3",
    size = 2.5,
    shape = 16,
    fill = "white",
    stroke = 1.2
  ) +
  geom_text(
    data = kaks_data,
    aes(x = Class, y = KaKs_mean_adjust, label = sprintf("%.2f", KaKs_mean)),
    color = "#1E395F",
    size = 4,
    vjust = -1.2,
    fontface = "bold"
  ) +
  
  # 双轴设置
  scale_y_continuous(
    name = "Ka",
    limits = main_y_range,
    expand = c(0, 0),
    sec.axis = sec_axis(~. / trans_coef, name = "Ks")
  ) +
  
  # 图例颜色配置（包含Ka、Ks、Ka/Ks）
  scale_fill_manual(
    values = c("Ka" = "#5B859E", "Ks" = "#1E395F"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Ka/Ks" = "#DF8D71"),
    name = NULL
  ) +
  
  # 主题（图例位置在右上角）
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 1, vjust = 1, size = 14),
    axis.title.y = element_text(size = 16, color = "black", margin = margin(r = 10)),
    
    axis.title.y.right = element_text(size = 16, color = "black", margin = margin(l = 10),angle = 90,vjust = 0.5),
    # 图例位置与样式（核心修改）
    legend.position = c(0.95, 0.95),  # 右上角坐标
    legend.justification = c(1, 1),   # 右上角对齐
    legend.margin = margin(2, 2, 2, 2),
    legend.key.size = unit(0.6, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.8)
    # 加粗边框线
        # 保持指定的边框颜色
    
  ) 
linesaa
ggsave("linesaa.tif",plot = linesaa,width = 2400,height = 1200,units = "px")

###save class ka/ks data
write.csv(kaks_data,file= "kaksdata.csv",quote =T)

#save Class ka ks data 
write.csv(stats_data,file= "stats_data.csv",quote =T)




