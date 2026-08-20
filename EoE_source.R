color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


plot_umap1 <- function(seurat_obj, 
                                 group.by = "celltype2", 
                                 reduction = "umap", 
                                 label = TRUE, 
                                 cols = NULL,
                                 arrow_angle = 20,
                                 arrow_length = 0.1,
                                 arrow_type = "closed",
                                 umap_xlim = c(-0.3, 1),
                                 umap_ylim = c(-0.3, 1),
                                 arrow_segment_size = 0.5,
                                 arrow_lineend = "round",
                                 # 修改参数名称为 arrow_text_size，表示文本大小
                                 arrow_text_size = 4,
                                 umap_coord_width = 0.2,
                                 umap_coord_height = 0.2,
                                 main_plot_scale = 0.9) {
  
  # 1. 构造箭头对象
  my_arrow <- arrow(angle = arrow_angle, 
                    type = arrow_type, 
                    length = unit(arrow_length, "npc"))
  
  # 2. 绘制 UMAP 坐标系箭头示意图
  umap_coord <- ggplot(tibble(
    group = c("UMAP1", "UMAP2"),
    x     = c(0, 0), 
    xend  = c(1, 0),
    y     = c(0, 0), 
    yend  = c(0, 1),
    lx    = c(0.5, -0.15), 
    ly    = c(-0.15, 0.5),
    angle = c(0, 90)
  )) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, group = group),
                 arrow = my_arrow, 
                 linewidth = arrow_segment_size,    # 使用 linewidth 替换 size
                 lineend = arrow_lineend) +
    geom_text(aes(x = lx, y = ly, label = group, angle = angle),
              size = arrow_text_size) +         # 使用 size 设置文本大小
    theme_void() +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim)
  
  color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
  
  
  # 3. 绘制主 UMAP 图，使用 Seurat 的 DimPlot
  p_dim <- DimPlot(object = seurat_obj, 
                   reduction = reduction, 
                   group.by = group.by,
                   label = label, 
                   cols = cols)
  
  # 4. 利用 cowplot 叠加两幅图
  combined_plot <- ggdraw() +
    draw_plot(
      p_dim + NoAxes(),   # 去除主图的坐标轴
      scale = main_plot_scale
    ) +
    draw_plot(
      umap_coord,
      x = 0, y = 0,
      width = umap_coord_width, height = umap_coord_height
    )
  
  return(combined_plot)
}


plot_umap2_exp2 <- function(seurat_obj, 
                           gene_name, 
                           min_value = -4, 
                           max_value = 4,
                           arrow_angle = 20,
                           arrow_length = 0.1,
                           arrow_segment_size = 0.5,
                           arrow_lineend = "round",   # 用于 geom_segment() 的 lineend
                           arrow_type = "closed",     # 用于 arrow() 的 type 参数，必须是 "open" 或 "closed"
                           arrow_text_size = 4,
                           umap_xlim = c(-0.3, 1),
                           umap_ylim = c(-0.3, 1),
                           umap_coord_width = 0.2,
                           umap_coord_height = 0.2,
                           title_size = 14) {
  # 加载必要的包
  library(Seurat)
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(tibble)
  library(scales)
  
  # -------------------------
  # 1. 数据准备与处理
  # -------------------------
  
  # 获取基因的 z-score 数据（scale.data 层）
  z_scores <- FetchData(seurat_obj, vars = gene_name, layer = "scale.data")
  
  # 获取 UMAP 坐标，并整理为 data.frame
  umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_data) <- c("umap_1", "umap_2")
  
  # 利用 FeaturePlot 获取初步数据
  p_tmp <- FeaturePlot(seurat_obj, features = gene_name, min.cutoff = min_value, max.cutoff = max_value)
  data <- p_tmp$data
  
  # 将 UMAP 坐标信息合并到 data 中
  data$umap_1 <- umap_data[rownames(data), "umap_1"]
  data$umap_2 <- umap_data[rownames(data), "umap_2"]
  
  # 手动修剪表达数据：将表达值限制在 [min_value, max_value] 内
  data[[gene_name]][data[[gene_name]] < min_value] <- min_value
  data[[gene_name]][data[[gene_name]] > max_value] <- max_value
  
  # 根据表达量排序（便于低表达的点在底层）
  data <- data[order(data[[gene_name]]), ]
  
  # 添加元数据：提取 cell 对应的 condition（这里存于 meta.data 中的 "condition" 字段）和 celltype 信息
  cell_ids <- rownames(data)
  data$Timepoint <- seurat_obj@meta.data[cell_ids, "condition"]
  data$celltype <- seurat_obj@meta.data[cell_ids, "celltype"]

  
  # -------------------------
  # 2. 绘制表达量散点图（去除坐标轴）
  # -------------------------
  
  # 先绘制背景黑色点（这里以 celltype 作为映射，仅用作背景展示黑点）
  p0 <- ggplot(data, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 1.5, color = "black")
  
  # 利用 ggnewscale 切换新的颜色映射
  p1 <- p0 + ggnewscale::new_scale_color()
  
  # 添加表达量图层；此处不添加标题，后续通过 cowplot 单独绘制标题
  p_final <- p1 +
    geom_point(aes(color = .data[[gene_name]]), size = 0.8, alpha = 1) +
    scale_color_gradientn(colors = c('#2E5A87', 'gray97', '#A90C38'),
                          limits = c(min_value, max_value),
                          oob = squish,
                          breaks = c(min_value, 0, max_value)) +
    theme_void() +
    theme(
      legend.position = "right"
    )
  
  # -------------------------
  # 3. 构造箭头示意图
  # -------------------------
  
  # 使用正确的 arrow_type 参数来构造箭头对象
  my_arrow <- arrow(angle = arrow_angle, 
                    type = arrow_type, 
                    length = unit(arrow_length, "npc"))
  
  umap_coord <- ggplot(tibble(
    group = c("UMAP1", "UMAP2"),
    x     = c(0, 0),
    xend  = c(1, 0),
    y     = c(0, 0),
    yend  = c(0, 1),
    lx    = c(0.5, -0.15),
    ly    = c(-0.15, 0.5),
    angle = c(0, 90)
  )) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, group = group),
                 arrow = my_arrow, 
                 linewidth = arrow_segment_size, 
                 lineend = arrow_lineend) +
    geom_text(aes(x = lx, y = ly, label = group, angle = angle), 
              size = arrow_text_size) +
    theme_void() +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim)
  
  # -------------------------
  # 4. 组合主图、箭头图以及标题
  # -------------------------
  
  # 使用 cowplot::ggdraw 组合各图层，并利用 draw_text 在左上角绘制标题
  combined_plot <- ggdraw() +
    # 主图全幅显示
    draw_plot(p_final, x = 0, y = 0, width = 1, height = 1) +
    # 叠加箭头示意图（默认放置于左下角，可根据需要调整位置）
    draw_plot(umap_coord, x = 0, y = 0, width = umap_coord_width, height = umap_coord_height) +
    # 添加标题：位于左上角，x 与 y 坐标采用相对比例（例如 x = 0.01, y = 0.99）
    draw_text(text = gene_name, x = 0.01, y = 0.99, hjust = 0, vjust = 1,
              size = title_size, fontface = "bold")
  
  return(combined_plot)
}


plot_umap2_exp <- function(seurat_obj, 
                           gene_name, 
                           arrow_angle = 20,
                           arrow_length = 0.1,
                           arrow_segment_size = 0.5,
                           arrow_lineend = "round",   # 用于 geom_segment() 的 lineend
                           arrow_type = "closed",     # 用于 arrow() 的 type 参数，必须是 "open" 或 "closed"
                           arrow_text_size = 4,
                           umap_xlim = c(-0.3, 1),
                           umap_ylim = c(-0.3, 1),
                           umap_coord_width = 0.2,
                           umap_coord_height = 0.2,
                           title_size = 14) {
  # 加载必要的包
  library(Seurat)
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(tibble)
  library(scales)
  
  # -------------------------
  # 1. 数据准备与处理
  # -------------------------
  
  # 获取 UMAP 坐标
  umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_data) <- c("umap_1", "umap_2")
  
  # 利用 FeaturePlot 获取基因表达数据
  p_tmp <- FeaturePlot(seurat_obj, features = gene_name)
  data <- p_tmp$data
  
  # 将 UMAP 坐标信息合并到 data 中
  data$umap_1 <- umap_data[rownames(data), "umap_1"]
  data$umap_2 <- umap_data[rownames(data), "umap_2"]
  
  # 获取表达数据的最大值（自动判断，并取1位小数）
  max_value <- max(data[[gene_name]], na.rm = TRUE)
  max_value <- round(max_value, 1)
  
  # 根据表达量排序（便于低表达的点在底层）
  data <- data[order(data[[gene_name]]), ]
  
  # 添加元数据（假设 meta.data 中包含 "condition" 和 "celltype"）
  cell_ids <- rownames(data)
  data$Timepoint <- seurat_obj@meta.data[cell_ids, "condition"]
  data$celltype <- seurat_obj@meta.data[cell_ids, "celltype"]
  
  # -------------------------
  # 2. 绘制表达量散点图（去除坐标轴）
  # -------------------------
  
  # 先绘制背景黑色点（以 celltype 作为映射，仅用作背景展示黑点）
  p0 <- ggplot(data, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 1.2, color = "black")
  
  # 利用 ggnewscale 切换新的颜色映射
  p1 <- p0 + ggnewscale::new_scale_color()
  
  # 添加表达量图层，并设置颜色映射：0 显示为 "0"，最大值保留1位小数
  p_final <- p1 +
    geom_point(aes(color = .data[[gene_name]]), size = 0.4, alpha = 1) +
    scale_color_gradientn(colors = c('gray97', '#A90C38'),
                          limits = c(0, max_value),
                          oob = squish,
                          breaks = c(0, max_value),
                          labels = c("0", sprintf("%.1f", max_value))) +
    theme_void() +
    theme(
      legend.position = "right"
    )
  
  # -------------------------
  # 3. 构造箭头示意图
  # -------------------------
  
  # 使用正确的 arrow_type 参数来构造箭头对象
  my_arrow <- arrow(angle = arrow_angle, 
                    type = arrow_type, 
                    length = unit(arrow_length, "npc"))
  
  umap_coord <- ggplot(tibble(
    group = c("UMAP1", "UMAP2"),
    x     = c(0, 0),
    xend  = c(1, 0),
    y     = c(0, 0),
    yend  = c(0, 1),
    lx    = c(0.5, -0.15),
    ly    = c(-0.15, 0.5),
    angle = c(0, 90)
  )) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, group = group),
                 arrow = my_arrow, 
                 linewidth = arrow_segment_size, 
                 lineend = arrow_lineend) +
    geom_text(aes(x = lx, y = ly, label = group, angle = angle), 
              size = arrow_text_size) +
    theme_void() +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim)
  
  # -------------------------
  # 4. 组合主图、箭头图以及标题
  # -------------------------
  
  # 使用 cowplot::ggdraw 组合各图层，并利用 draw_text 在左上角绘制标题
  combined_plot <- ggdraw() +
    # 主图全幅显示
    draw_plot(p_final, x = 0, y = 0, width = 1, height = 1) +
    # 叠加箭头示意图（默认放置于左下角，可根据需要调整位置）
    draw_plot(umap_coord, x = -0.06, y = 0, width = umap_coord_width, height = umap_coord_height) +
    # 添加标题：位于左上角，x 与 y 坐标采用相对比例（例如 x = 0.01, y = 0.99）
    draw_text(text = gene_name, x = 0.01, y = 0.99, hjust = 0, vjust = 1,
              size = title_size, fontface = "bold")
  
  return(combined_plot)
}



plot_umap4_facet <- function(seurat_obj, 
                             gene_name, 
                             arrow_angle = 20,
                             arrow_length = 0.1,
                             arrow_segment_size = 0.5,
                             arrow_lineend = "round",   # 用于 geom_segment() 的 lineend
                             arrow_type = "closed",     # 用于 arrow() 的 type 参数，必须是 "open" 或 "closed"
                             arrow_text_size = 4,
                             umap_xlim = c(-0.3, 1),
                             umap_ylim = c(-0.3, 1),
                             umap_coord_width = 0.2,
                             umap_coord_height = 0.2,
                             title_size = 14) {
  # 加载必要的包
  library(Seurat)
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(tibble)
  library(scales)
  
  # -------------------------
  # 1. 数据准备与处理
  # -------------------------
  
  # 获取 UMAP 坐标
  umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_data) <- c("umap_1", "umap_2")
  
  # 利用 FeaturePlot 获取基因表达数据
  p_tmp <- FeaturePlot(seurat_obj, features = gene_name)
  data <- p_tmp$data
  
  # 将 UMAP 坐标信息合并到 data 中
  data$umap_1 <- umap_data[rownames(data), "umap_1"]
  data$umap_2 <- umap_data[rownames(data), "umap_2"]
  
  # 获取表达数据的最大值（自动判断，并取1位小数）
  max_value <- max(data[[gene_name]], na.rm = TRUE)
  max_value <- round(max_value, 1)
  
  # 根据表达量排序（便于低表达的点在底层）
  data <- data[order(data[[gene_name]]), ]
  
  # 添加元数据（假设 meta.data 中包含 "condition" 和 "celltype"）
  cell_ids <- rownames(data)
  data$Timepoint <- seurat_obj@meta.data[cell_ids, "condition"]
  data$celltype <- seurat_obj@meta.data[cell_ids, "celltype"]
  
  # 确保 Timepoint 顺序为 Health 在左, EoE 在右
  data$Timepoint <- factor(data$Timepoint, levels = c("Health", "EoE"))
  
  # -------------------------
  # 2. 绘制表达量散点图（去除坐标轴）
  # -------------------------
  
  # 先绘制背景黑色点（这里以 celltype 作为映射，仅用作背景展示黑点）
  p0 <- ggplot(data, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 1.2, color = "black")
  
  # 利用 ggnewscale 切换新的颜色映射
  p1 <- p0 + ggnewscale::new_scale_color()
  
  # 添加表达量图层，并采用 facet_wrap 分面，同时为 facet strip 添加灰色底色和黑色边框
  p_final <- p1 +
    geom_point(aes(color = .data[[gene_name]]), size = 0.4, alpha = 1) +
    scale_color_gradientn(colors = c('gray97', '#A90C38'),
                          limits = c(0, max_value),
                          oob = squish,
                          breaks = c(0, max_value),
                          labels = c("0", sprintf("%.1f", max_value))) +
    facet_wrap(~ Timepoint, ncol = 2) +  # 这里保证 facet 分面顺序正确
    theme_void() +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(color = "black")
    )
  
  # -------------------------
  # 3. 构造箭头示意图
  # -------------------------
  
  # 使用正确的 arrow_type 参数来构造箭头对象
  my_arrow <- arrow(angle = arrow_angle, 
                    type = arrow_type, 
                    length = unit(arrow_length, "npc"))
  
  umap_coord <- ggplot(tibble(
    group = c("UMAP1", "UMAP2"),
    x     = c(0, 0),
    xend  = c(1, 0),
    y     = c(0, 0),
    yend  = c(0, 1),
    lx    = c(0.5, -0.15),
    ly    = c(-0.15, 0.5),
    angle = c(0, 90)
  )) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, group = group),
                 arrow = my_arrow, 
                 linewidth = arrow_segment_size, 
                 lineend = arrow_lineend) +
    geom_text(aes(x = lx, y = ly, label = group, angle = angle), 
              size = arrow_text_size) +
    theme_void() +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim)
  
  # -------------------------
  # 4. 组合主图、箭头图以及标题
  # -------------------------
  
  # 使用 cowplot::ggdraw 组合各图层，并利用 draw_text 在左上角绘制标题
  combined_plot <- ggdraw() +
    # 主图全幅显示
    draw_plot(p_final, x = 0, y = 0, width = 1, height = 1) +
    # 叠加箭头示意图（默认放置于左下角，可根据需要调整位置）
    draw_plot(umap_coord, x = -0.065, y = -0.07, width = umap_coord_width, height = umap_coord_height) +
    # 添加标题：位于左上角，x 与 y 坐标采用相对比例（例如 x = 0.01, y = 0.99）
    draw_text(text = gene_name, x = 0.01, y = 0.93, hjust = 0, vjust = 1,
              size = title_size, fontface = "bold")
  
  return(combined_plot)
}





plot_umap4_facet2 <- function(seurat_obj, 
                             gene_name, 
                             arrow_angle = 20,
                             arrow_length = 0.1,
                             arrow_segment_size = 0.5,
                             arrow_lineend = "round",   # 用于 geom_segment() 的 lineend
                             arrow_type = "closed",     # 用于 arrow() 的 type 参数，必须是 "open" 或 "closed"
                             arrow_text_size = 4,
                             umap_xlim = c(-0.3, 1),
                             umap_ylim = c(-0.3, 1),
                             umap_coord_width = 0.2,
                             umap_coord_height = 0.2,
                             title_size = 14) {
  # 加载必要的包
  library(Seurat)
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(tibble)
  library(scales)
  
  # -------------------------
  # 1. 数据准备与处理
  # -------------------------
  
  # 获取 UMAP 坐标
  umap_data <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_data) <- c("umap_1", "umap_2")
  
  # 利用 FeaturePlot 获取基因表达数据
  p_tmp <- FeaturePlot(seurat_obj, features = gene_name)
  data <- p_tmp$data
  
  # 将 UMAP 坐标信息合并到 data 中
  data$umap_1 <- umap_data[rownames(data), "umap_1"]
  data$umap_2 <- umap_data[rownames(data), "umap_2"]
  
  # 获取表达数据的最大值（自动判断，并取1位小数）
  max_value <- max(data[[gene_name]], na.rm = TRUE)
  max_value <- round(max_value, 1)
  
  # 根据表达量排序（便于低表达的点在底层）
  data <- data[order(data[[gene_name]]), ]
  
  # 添加元数据（假设 meta.data 中包含 "condition" 和 "celltype"）
  cell_ids <- rownames(data)
  data$Timepoint <- seurat_obj@meta.data[cell_ids, "condition"]
  data$celltype <- seurat_obj@meta.data[cell_ids, "celltype"]
  
  # 确保 Timepoint 顺序为 Health 在左, EoE 在右
  data$Timepoint <- factor(data$Timepoint, levels = c("Health", "EoE"))
  
  # -------------------------
  # 2. 绘制表达量散点图（去除坐标轴）
  # -------------------------
  
  # 先绘制背景黑色点（这里以 celltype 作为映射，仅用作背景展示黑点）
  p0 <- ggplot(data, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 0, color = "black")
  
  # 利用 ggnewscale 切换新的颜色映射
  p1 <- p0 + ggnewscale::new_scale_color()
  
  # 添加表达量图层，并采用 facet_wrap 分面，同时为 facet strip 添加灰色底色和黑色边框
  p_final <- p1 +
    geom_point(aes(color = .data[[gene_name]]), size = 0.4, alpha = 1) +
    scale_color_gradientn(colors = c('gray97', '#A90C38'),
                          limits = c(0, max_value),
                          oob = squish,
                          breaks = c(0, max_value),
                          labels = c("0", sprintf("%.1f", max_value))) +
    facet_wrap(~ Timepoint, ncol = 2) +  # 这里保证 facet 分面顺序正确
    theme_void() +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(color = "black")
    )
  
  # -------------------------
  # 3. 构造箭头示意图
  # -------------------------
  
  # 使用正确的 arrow_type 参数来构造箭头对象
  my_arrow <- arrow(angle = arrow_angle, 
                    type = arrow_type, 
                    length = unit(arrow_length, "npc"))
  
  umap_coord <- ggplot(tibble(
    group = c("UMAP1", "UMAP2"),
    x     = c(0, 0),
    xend  = c(1, 0),
    y     = c(0, 0),
    yend  = c(0, 1),
    lx    = c(0.5, -0.15),
    ly    = c(-0.15, 0.5),
    angle = c(0, 90)
  )) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, group = group),
                 arrow = my_arrow, 
                 linewidth = arrow_segment_size, 
                 lineend = arrow_lineend) +
    geom_text(aes(x = lx, y = ly, label = group, angle = angle), 
              size = arrow_text_size) +
    theme_void() +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim)
  
  # -------------------------
  # 4. 组合主图、箭头图以及标题
  # -------------------------
  
  # 使用 cowplot::ggdraw 组合各图层，并利用 draw_text 在左上角绘制标题
  combined_plot <- ggdraw() +
    # 主图全幅显示
    draw_plot(p_final, x = 0, y = 0, width = 1, height = 1) +
    # 叠加箭头示意图（默认放置于左下角，可根据需要调整位置）
    draw_plot(umap_coord, x = -0.065, y = -0.07, width = umap_coord_width, height = umap_coord_height) +
    # 添加标题：位于左上角，x 与 y 坐标采用相对比例（例如 x = 0.01, y = 0.99）
    draw_text(text = gene_name, x = 0.01, y = 0.93, hjust = 0, vjust = 1,
              size = title_size, fontface = "bold")
  
  return(combined_plot)
}






plot_violin_1vs1 <- function(test.seu, gene, subset = NULL) {
  # 如果传入了subset，则根据celltype2进行子集筛选
  if (!is.null(subset)) {
    test.seu <- subset(test.seu, celltype2 == subset)
  }
  
  # 如果条件中包含"Health"和"EoE"，重新设置因子水平顺序，使"Health"显示在左边
  if(all(c("Health", "EoE") %in% unique(test.seu$condition))) {
    test.seu$condition <- factor(test.seu$condition, levels = c("Health", "EoE"))
  }
  
  # 确保感兴趣的基因在数据中存在
  genes_of_interest <- intersect(gene, rownames(test.seu))
  if (length(genes_of_interest) == 0) {
    stop(paste("基因", gene, "不在对象中！"))
  }
  
  # 提取该基因的表达数据，并添加condition分组信息
  expression_data <- FetchData(test.seu, vars = genes_of_interest)
  expression_data$condition <- test.seu$condition
  
  # 配置颜色（此处生成的grad_ungroup未在后续绘图中使用，可根据需要调整）
  colours <- scales::viridis_pal()(10)
  grad_ungroup <- linearGradient(colours, group = FALSE)
  
  # 创建小提琴图，按condition分组
  p1 <- VlnPlot(
    object = test.seu,
    features = gene,
    group.by = "condition",
    pt.size = 0,
    assay = "RNA",
    cols = c("Health" = "#b7b7b7", "EoE" = "#B72230")
  ) +
    theme_prism(
      palette = "colors",
      base_fontface = "plain",
      base_family = "sans",
      base_size = 16,
      base_line_size = 0.8,
      axis_text_angle = 45
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste(gene, "Expression")) +
    ylim(-0.01, max(expression_data[[gene]], na.rm = TRUE) * 1.3)
  
  # 添加显著性检验（t检验），比较"Health"和"EoE"两组
  p2 <- p1 + geom_signif(
    comparisons = list(c("Health", "EoE")),
    map_signif_level = TRUE,  # 显示显著性水平
    test = "t.test",          # 使用t检验
    test.args = list(alternative = "two.sided", var.equal = FALSE),
    y_position = max(FetchData(test.seu, gene), na.rm = TRUE) * 1.1,
    tip_length = 0,
    textsize = 6,
    size = 0.6,
    color = "black"
  )
  
  # 如果输入了subset，则标题修改为 "gene in subset"
  if (!is.null(subset)) {
    p2 <- p2 + 
      ggtitle(paste(gene, "in", subset)) +
      theme(plot.title = element_text(size = 14))
  }
  
  
  # 计算每组的平均表达值
  mean_health <- mean(expression_data[expression_data$condition == "Health", gene], na.rm = TRUE)
  mean_eoe    <- mean(expression_data[expression_data$condition == "EoE", gene], na.rm = TRUE)
  
  # 修改横坐标标签：在Health和EoE后另起一行，显示平均值（保留一位小数）
  p2 <- p2 + scale_x_discrete(labels = c(
    "Health" = paste0("Health\n(", round(mean_health, 2), ")"),
    "EoE"    = paste0("EoE\n(", round(mean_eoe, 2), ")")
  ))
  
  # 返回最终图形
  return(p2)
}



plot_EoE <- function(seurat_obj,                # Seurat object containing the data
                     gene,                      # Gene name to be visualized
                     group.by = "celltype2",    # Metadata column used to group cells (default "celltype2")
                     subset = "Quiescent Basal",# Specific subset of cells to focus on (default "Quiescent Basal")
                     output_file = NULL,        # File name to save the output plot (default NULL, auto-generated)
                     output_format = "pdf",     # Output format: "pdf" or "png" (default "pdf")
                     width = 50,                # Width of the saved plot (default 50)
                     height = 30,               # Height of the saved plot (default 30)
                     units = "cm") {            # Units for width and height (default "cm")
  
  # Convert output_format to lowercase to ensure consistency
  output_format <- tolower(output_format)  # Ensure format is in lowercase
  
  # Check if output_format is either "pdf" or "png" and stop if not
  if (!(output_format %in% c("pdf", "png"))) {  # Validate output format
    stop("output_format MUST BE 'pdf' or 'png'")  # Stop execution if invalid format is provided
  }
  
  # If no output_file is provided, generate a default file name based on gene and subset
  if (is.null(output_file)) {  # Check if output_file is NULL
    if (!is.null(subset) && subset != "") {  # If subset is provided and not an empty string
      output_file <- paste0("[EoE]", gene, "_exp_{", subset, "}.", output_format)  # Construct file name including subset
    } else {
      output_file <- paste0("[EoE]", gene, "_exp.", output_format)  # Construct file name without subset
    }
  }
  
  # Source external R file containing custom plotting functions; ensure EoE_source.R is in the working directory
  source('EoE_source.R')  # Load additional plotting functions from external file
  
  # Define a color scheme vector for the plots
  color_scheme <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")  # Set up a vector of colors
  
  # Define an alternative color vector (can be used for specific cell types)
  color_cc <- c("#E15759", "#F0C420", "#1965B0")  # Alternative color scheme
  
  # 1. Create a UMAP clustering plot using the specified grouping variable (group.by)
  u1 <- plot_umap1(seurat_obj,                  # Pass the Seurat object
                   group.by = group.by,            # Group cells by this metadata field
                   reduction = "umap",             # Use the UMAP reduction
                   label = TRUE,                   # Add labels to clusters
                   cols = color_scheme)            # Apply the defined color scheme
  
  # 2. Create a UMAP plot showing the expression of the specified gene
  u2 <- plot_umap2_exp(seurat_obj = seurat_obj,  # Pass the Seurat object
                       gene_name = gene)           # Specify the gene to visualize
  
  # 3. Create a gene expression density plot on the UMAP embedding
  u3 <- plot_density(seurat_obj,                # Seurat object
                     gene,                      # Gene to plot
                     size = 0.8,                # Point size for density estimation
                     joint = FALSE,             # Not plotting joint density
                     reduction = "umap",        # Use UMAP coordinates
                     combine = TRUE,            # Combine layers if applicable
                     pal = "inferno")  +        # Use the "inferno" color palette
    theme_void()                              # Remove background and axes for a clean look
  
  # 4. Create a faceted UMAP plot displaying gene expression across facets
  u4 <- plot_umap4_facet(seurat_obj = seurat_obj,  # Seurat object
                         gene_name = gene)           # Gene to display
  
  # 5. Generate a violin plot of gene expression split by condition (e.g., Health vs EoE) and grouped by group.by
  v2 <- VlnPlot(seurat_obj,                     # Seurat object
                features = gene,                # Gene to plot
                stack = FALSE,                  # Do not stack violins
                pt.size = 0,                    # Hide individual data points
                flip = FALSE,                   # Do not flip coordinates
                add.noise = TRUE,               # Add noise to the data points (jitter)
                split.by = 'condition',         # Split violins by condition
                group.by = group.by,            # Group violins by specified metadata field
                cols = c("#C73E3A", "#b7b7b7"),  # Colors for the two conditions
                split.plot = TRUE)              # Create a split plot
  
  # Adjust the order of the split factor so that "Health" comes before "EoE"
  v2$data$split <- factor(v2$data$split, levels = c("Health", "EoE"))  # Reorder factor levels for split variable
  
  # Add a border to each panel of the violin plot
  v2 <- v2 + theme(panel.border = element_rect(color = "black", 
                                               fill = NA, 
                                               size = 1))  # Apply a black border to the panel
  
  # 6. Create a violin plot comparing gene expression across all cells
  v3 <- plot_violin_1vs1(seurat_obj, gene = gene)  # Violin plot for all cells
  
  # 7. Create a violin plot for a specified subset of cells; if subset is NULL, use the full set
  if (!is.null(subset)) {  # Check if a subset is provided
    v4 <- plot_violin_1vs1(seurat_obj, gene = gene, subset = subset)  # Create violin plot for the specified subset
  } else {
    v4 <- plot_violin_1vs1(seurat_obj, gene = gene)  # Otherwise, create violin plot for all cells
  }
  
  # Overlay a crossbar on the subset violin plot to indicate the mean expression value (displayed as a dotted line)
  v4 <- v4 + stat_summary(fun = "mean",        # Calculate the mean of the data
                          geom = "crossbar",     # Display it as a crossbar
                          width = 0.6,           # Set the width of the crossbar
                          colour = "black",      # Set the crossbar color to black
                          size = 0.2,            # Set the line thickness
                          linetype = "dotted",   # Use a dotted line style
                          show.legend = FALSE)   # Do not show a legend for the crossbar
  
  # Define the layout for combining plots using patchwork; each letter represents a plot
  design <- "
AAAAAAAAAAAAAABBBBBBBBBCCCCCCC
AAAAAAAAAAAAAABBBBBBBBBCCCCCCC
AAAAAAAAAAAAAABBBBBBBBBCCCCCCC
AAAAAAAAAAAAAABBBBBBBBBCCCCCCC
DDDDDDDDDDDDDDDDDDDDDDEEEEEEEE
DDDDDDDDDDDDDDDDDDDDDDFFFFGGGG
"  # Layout design string specifying where each plot goes
  
  # Combine all plots into one layout using patchwork and add panel tags
  p4 <- u1 + u2 + u3 + u4 + v2 + v3 + v4 +         # Combine the individual plots
    plot_layout(design = design) +                  # Apply the custom layout design
    plot_annotation(tag_levels = "A")               # Annotate each panel with letters (A, B, C, …)
  
  # Print the combined plot to the active graphic device
  print(p4)  # Output the final combined plot
  
  # Save the combined plot to a file using the specified output format and dimensions
  ggsave(filename = output_file,                     # File name for saving the plot
         plot = p4,                                 # Plot object to be saved
         device = output_format,                    # Device to use based on output_format ("pdf" or "png")
         dpi = 'retina',                            # Set high DPI for Retina displays
         width = width,                             # Width of the saved plot
         height = height,                           # Height of the saved plot
         units = units)                             # Units for width and height
  
  # Return the combined plot object so that it can be further manipulated if needed
  return(p4)  # Return the patchwork plot object
}








# neb_custom.R
# ------------

# 0) 加载依赖
library(Nebulosa)
library(ggplot2)
library(ggrastr)
library(patchwork)

# 1) 导入 Nebulosa 内部辅助函数
.validate_dimensions    <- Nebulosa:::.validate_dimensions
.search_dimensions      <- Nebulosa:::.search_dimensions
.extract_feature_data   <- Nebulosa:::.extract_feature_data
calculate_density       <- Nebulosa:::calculate_density

# 2) 重写泛型，让 pal 支持更多选项
setGeneric("plot_density", function(object, features, slot = NULL,
                                    joint = FALSE, reduction = NULL,
                                    dims = c(1, 2),
                                    method = c("ks", "wkde"),
                                    adjust = 1, size = 1, shape = 16,
                                    combine = TRUE,
                                    pal = c(
                                      "viridis","magma","cividis","inferno","plasma",
                                      "rocket","mako","turbo",     
                                      "A","B","C","D","E","F","G","H"  
                                    ),
                                    raster = TRUE,
                                    ...) standardGeneric("plot_density"))

# 3) 重写 Seurat 方法
setMethod("plot_density", signature("Seurat"),
          function(object, features, slot = NULL, joint = FALSE, reduction = NULL,
                   dims = c(1, 2), method = c("ks", "wkde"), adjust = 1,
                   size = 1, shape = 16, combine = TRUE,
                   pal = c(
                     "viridis","magma","cividis","inferno","plasma",
                     "rocket","mako","turbo","A","B","C","D","E","F","G","H"
                   ),
                   raster = TRUE, ...) {
            
            .validate_dimensions(dims)
            if (!is.null(reduction) && !reduction %in% Reductions(object)) {
              stop("No reduction named '", reduction, "' found in object")
            }
            reduction_list <- Reductions(object)
            if (length(reduction_list) == 0) stop("No reduction has been computed!")
            if (is.null(reduction)) reduction <- reduction_list[length(reduction_list)]
            cell_embeddings <- as.data.frame(Embeddings(object[[reduction]]))
            cell_embeddings <- .search_dimensions(dims, cell_embeddings, reduction)
            
            if (is.null(slot)) slot <- "data"
            method <- match.arg(method)
            metadata <- object[[]]
            if (all(features %in% colnames(metadata))) {
              vars <- metadata[, features, drop = FALSE]
            } else {
              exp_data <- FetchData(object, vars = features, slot = slot)
              vars <- .extract_feature_data(exp_data, features)
            }
            
            .plot_final_density(vars, cell_embeddings, features,
                                joint, method, adjust,
                                shape, size, pal, combine, raster, ...)
          })

# 4) 重写组合逻辑
.plot_final_density <- function(vars, cell_embeddings, features,
                                joint, method, adjust,
                                shape, size, pal, combine, raster, ...) {
  dim_names <- colnames(cell_embeddings)
  
  if (ncol(vars) > 1) {
    res <- apply(vars, 2, calculate_density, cell_embeddings, method, adjust)
    p_list <- mapply(plot_density_, as.list(as.data.frame(res)), colnames(res),
                     MoreArgs = list(cell_embeddings, dim_names,
                                     shape, size, "Density", pal, raster, ...),
                     SIMPLIFY = FALSE)
    if (combine) {
      patchwork::wrap_plots(p_list, byrow = TRUE)
    } else {
      p_list
    }
  } else {
    z <- calculate_density(vars[,1], cell_embeddings, method, adjust)
    plot_density_(z, features, cell_embeddings, dim_names,
                  shape, size, "Density", pal, raster, ...)
  }
}

# 5) 重写底层绘图函数，支持单色和多色
plot_density_ <- function(z, feature, cell_embeddings,
                          dim_names, shape, size, legend_title,
                          pal = c(
                            "viridis","magma","cividis","inferno","plasma",
                            "rocket","mako","turbo","A","B","C","D","E","F","G","H"
                          ),
                          raster, ...) {
  p <- ggplot2::ggplot(
    data.frame(cell_embeddings, feature = z),
    ggplot2::aes_string(dim_names[1], dim_names[2], color = "feature")
  ) +
    ggplot2::geom_point(shape = shape, size = size) +
    ggplot2::xlab(gsub("_", " ", dim_names[1])) +
    ggplot2::ylab(gsub("_", " ", dim_names[2])) +
    ggplot2::ggtitle(feature) +
    ggplot2::labs(color = legend_title) +
    ggplot2::theme(
      text             = ggplot2::element_text(size = 14),
      panel.background = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(color = "black"),
      axis.text.y      = ggplot2::element_text(color = "black"),
      axis.line        = ggplot2::element_line(size = 0.25),
      strip.background = ggplot2::element_rect(color = "black", fill = "#ffe5cc")
    )
  
  pal <- match.arg(pal)
  # 调用 ggplot2 中的 scale_colour_viridis_c
  if (length(pal) == 1 && pal %in% c("A","B","C","D","E","F","G","H",
                                     "viridis","magma","cividis",
                                     "inferno","plasma","rocket","mako","turbo")) {
    p2 <- p + ggplot2::scale_colour_viridis_c(option = pal, ...)
  } else {
    p2 <- p + ggplot2::scale_color_gradientn(colors = pal, ...)
  }
  
  if (raster) {
    ggrastr::rasterise(p2, dpi = 300)
  } else {
    p2
  }
}




plot_umap_density_ccm <- function(seurat_obj,
                                  gene,
                                  reduction    = "umap",
                                  dims         = c(1,2),
                                  method       = "ks",      # 默认改为 ks
                                  adjust       = 1,
                                  shape        = 16,
                                  size         = 0.8,
                                  color_ccm    = c('gray85','#b40001'),
                                  raster       = TRUE,
                                  legend_title = "Density") {
  # 1) 提取 UMAP 坐标
  if (! reduction %in% Seurat::Reductions(seurat_obj)) {
    stop("Reduction '", reduction, "' not found in object.")
  }
  coords <- Seurat::Embeddings(seurat_obj[[reduction]])[, dims]
  colnames(coords) <- c("UMAP_1","UMAP_2")
  
  # 2) 提取基因表达值
  if (gene %in% colnames(seurat_obj@meta.data)) {
    expr <- seurat_obj@meta.data[[gene]]
  } else {
    expr <- Seurat::FetchData(seurat_obj, vars = gene, slot = "data")[[gene]]
  }
  
  # 3) 计算核密度 (ks 方法)
  z <- Nebulosa:::calculate_density(expr, coords, method, adjust)
  
  # 4) 组装并排序：先画低密度，后画高密度
  df <- data.frame(UMAP_1 = coords[,"UMAP_1"],
                   UMAP_2 = coords[,"UMAP_2"],
                   density = z)
  df <- df[order(df$density), ]
  
  # 5) 绘图
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = density)) +
    {
      if (raster) {
        ggrastr::rasterise(geom_point(shape = shape, size = size), dpi = 300)
      } else {
        geom_point(shape = shape, size = size)
      }
    } +
    scale_color_gradientn(colors = color_ccm, name = legend_title) +
    theme_void() +
    ggtitle(paste0(gene, " density (", method, ")"))
  
  return(p)
}






plot_and_save_umap_density <- function(seurat_obj,
                                       gene,
                                       output_dir    = ".",
                                       prefix        = "umap_density",
                                       reduction     = "umap",
                                       method        = "ks",
                                       adjust        = 1,
                                       shape         = 16,
                                       size          = 0.8,
                                       color_ccm     = c('gray75','#b40001'),
                                       raster        = TRUE,
                                       legend_title  = "Density",
                                       width_cm      = 10,
                                       height_cm     = 10,
                                       dpi           = "retina",
                                       units         = "cm") {
  # 1) 生成密度图
  p <- plot_umap_density_ccm(
    seurat_obj = seurat_obj,
    gene       = gene,
    reduction  = reduction,
    method     = method,
    adjust     = adjust,
    shape      = shape,
    size       = size,
    color_ccm  = color_ccm,
    raster     = raster,
    legend_title = legend_title
  )
  
  # 2) 打印预览
  print(p)
  
  # 3) 构建文件名并保存
  fname <- file.path(
    output_dir,
    paste0(prefix, "_", gene, "_", method, ".pdf")
  )
  ggsave(
    filename = fname,
    plot     = p,
    device   = "pdf",
    dpi      = dpi,
    width    = width_cm,
    height   = height_cm,
    units    = units
  )
  
  message("Saved plot to: ", fname)
}






plot_and_save_density_umap <- function(seurat_obj,
                                       gene_name,
                                       output_prefix = "[EoE](den)",
                                       output_dir    = ".",
                                       pal           = "inferno",
                                       size          = 0.8,
                                       reduction     = "umap",
                                       joint         = FALSE,
                                       combine       = TRUE,
                                       width_cm      = 10,
                                       height_cm     = 10,
                                       dpi           = "retina",
                                       units         = "cm") {
  # 1) 生成密度图
  p <- plot_density(
    seurat_obj,
    gene_name,
    size      = size,
    joint     = joint,
    reduction = reduction,
    combine   = combine,
    pal       = pal
  ) + theme_void()
  
  # 2) 打印预览（可省略）
  print(p)
  
  # 3) 构建文件名并保存：在基因名后面加上下调色板的名称
  fname <- file.path(
    output_dir,
    paste0(output_prefix, gene_name, "_", pal, ".pdf")
  )
  ggsave(
    filename = fname,
    plot     = p,
    device   = "pdf",
    dpi      = dpi,
    width    = width_cm,
    height   = height_cm,
    units    = units
  )
  
  message("Saved: ", fname)
}







plot_umap4_facet_GSVA <- function(seurat_obj, 
                             gene_name, 
                             arrow_angle = 20,
                             arrow_length = 0.1,
                             arrow_segment_size = 0.5,
                             arrow_lineend = "round",
                             arrow_type = "closed",
                             arrow_text_size = 4,
                             umap_xlim = c(-0.3, 1),
                             umap_ylim = c(-0.3, 1),
                             umap_coord_width = 0.2,
                             umap_coord_height = 0.2,
                             title_size = 14) {
  # 必要包
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(scales)
  library(tibble)
  library(Seurat)
  library(grid)  # unit()
  
  # 1) UMAP 坐标
  umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_df) <- c("umap_1","umap_2")
  umap_df$cell <- rownames(umap_df)
  
  # 2) 提取表达/GSVA 分数
  if (gene_name %in% colnames(seurat_obj@meta.data)) {
    vals <- seurat_obj@meta.data[[gene_name]]
  } else if ("GSVA" %in% names(seurat_obj@assays) &&
             gene_name %in% rownames(GetAssayData(seurat_obj, assay="GSVA", layer="data"))) {
    mat  <- GetAssayData(seurat_obj, assay="GSVA", layer="data")
    vals <- mat[gene_name, ]
  } else {
    tmp  <- FetchData(seurat_obj, vars=gene_name, layer="data")
    vals <- tmp[[gene_name]]
  }
  exp_df <- data.frame(cell = names(vals), expr = vals, stringsAsFactors = FALSE)
  
  # 3) 合并
  df <- merge(umap_df, exp_df, by="cell", sort=FALSE)
  
  # 4) 分面信息
  df$Timepoint <- factor(
    seurat_obj@meta.data[df$cell, "condition"],
    levels = c("Health","EoE")
  )
  df$celltype <- seurat_obj@meta.data[df$cell, "celltype"]
  
  # 5) 排序
  df <- df[order(df$expr), ]
  
  # 6) 主图：背景黑点 + 表达点
  p0 <- ggplot(df, aes(x=umap_1, y=umap_2)) +
    geom_point(size=1.2, color="black")
  
  p_main <- p0 +
    new_scale_color() +
    geom_point(aes(color=expr), size=0.4) +
    scale_color_gradientn(
      colors = c('gray97','#A90C38'),
      limits = c(0, round(max(df$expr, na.rm=TRUE),1)),
      oob    = squish,
      breaks = c(0, round(max(df$expr, na.rm=TRUE),1)),
      labels = c("0", sprintf("%.1f", round(max(df$expr, na.rm=TRUE),1)))
    ) +
    facet_wrap(~ Timepoint, ncol=2) +
    theme_void() +
    theme(
      legend.position  = "right",
      strip.background = element_rect(fill="grey90", color="black"),
      strip.text       = element_text(color="black")
    )
  
  # 7) UMAP 坐标箭头 (显式映射，避免缺失)
  my_arrow <- arrow(angle=arrow_angle, type=arrow_type, length=unit(arrow_length,"npc"))
  coord_df <- tibble(
    group = c("UMAP1","UMAP2"),
    x     = c(0,0),  xend = c(1,0),
    y     = c(0,0),  yend = c(0,1),
    lx    = c(0.5,-0.15), ly = c(-0.15,0.5),
    angle = c(0,90)
  )
  p_coords <- ggplot(coord_df) +
    geom_segment(
      mapping   = aes(x=x, y=y, xend=xend, yend=yend, group=group),
      arrow     = my_arrow,
      linewidth = arrow_segment_size,
      lineend   = arrow_lineend
    ) +
    geom_text(
      mapping = aes(x=lx, y=ly, label=group, angle=angle),
      size    = arrow_text_size
    ) +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim) +
    theme_void()
  
  # 8) 合并并标题
  combined <- ggdraw() +
    draw_plot(p_main,   x=0, y=0, width=1, height=1) +
    draw_plot(p_coords, x=-0.065, y=-0.07,
              width=umap_coord_width, height=umap_coord_height) +
    draw_text(
      text      = gene_name,
      x         = 0.01, y = 0.93,
      hjust     = 0, vjust = 1,
      size      = title_size,
      fontface  = "bold"
    )
  
  return(combined)
}



plot_and_save_umap <- function(gene_name, seurat_obj, output_prefix = "[EoE]umap", output_dir = ".") {
  # 生成UMAP图
  umap_plot <- plot_umap2_exp(seurat_obj = seurat_obj, gene_name = gene_name)
  
  # 显示图像（可选）
  print(umap_plot)
  
  # 构建输出文件路径
  filename <- file.path(output_dir, paste0(output_prefix, "_", gene_name, ".pdf"))
  
  # 保存图像
  ggsave(filename = filename, plot = umap_plot, 
         device = "pdf", dpi = "retina", width = 10, height = 10, units = "cm")
}








plot_violin_by_patient <- function(seu_obj, gene, celltype_name, file_name) {
  # 1. Subset to the chosen celltype2 and assign to test.seu6 in global env
  test.seu6 <<- subset(seu_obj, subset = celltype2 == celltype_name)
  
  # 2. Compute patient ordering
  patient_df <- test.seu6@meta.data %>%
    distinct(patient, condition) %>%
    arrange(condition, patient)
  patient_order <- patient_df %>% pull(patient)
  test.seu6@meta.data$patient <- factor(test.seu6@meta.data$patient,
                                        levels = patient_order)
  
  # 3. Ensure RNA assay & fetch expression for ymax
  DefaultAssay(test.seu6) <- "RNA"
  expr_df  <- FetchData(test.seu6,
                        vars  = gene,
                        assay = "RNA",
                        slot  = "data")
  max_val  <- max(expr_df[[gene]], na.rm = TRUE)
  
  # 4. Build color vector based on number of Health/EoE patients
  nH <- sum(patient_df$condition == "Health")
  nE <- sum(patient_df$condition == "EoE")
  color_patient <- c(rep("#b7b7b7", nH), rep("#C73E3A", nE))
  
  # 5. Create the violin plot
  p <- VlnPlot(
    object    = test.seu6,
    features  = gene,
    group.by  = "patient",
    pt.size   = 0,
    assay     = "RNA",
    cols      = color_patient
  ) +
    theme_prism(
      palette         = "colors",
      base_fontface   = "plain",
      base_family     = "sans",
      base_size       = 16,
      base_line_size  = 0.8,
      axis_text_angle = 45
    ) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1),
      axis.title.x   = element_blank(),
      legend.position = "none"
    ) +
    labs(y = paste0(gene, " expression")) +
    coord_cartesian(ylim = c(0, max_val * 1.3))
  
  # 6. Save to PDF and return the plot
  ggsave(filename = file_name, plot = p,
         device = "pdf", dpi = 300,
         width = 20, height = 8, units = "cm")
  return(p)
}



squash_axis <- function(from, to, factor) { 
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}



squash_axis_symmetric2 <- function(from, to, factor) {
  trans <- function(x) {
    out <- x
    # 只对非 NA 且落在 (from, to) 以及 [to, +∞) 区间的点进行处理
    idx1 <- which(!is.na(x) & abs(x) > from & abs(x) < to)
    idx2 <- which(!is.na(x) & abs(x) >= to)
    # 压缩段
    out[idx1] <- sign(x[idx1]) * (from + (abs(x[idx1]) - from) / factor)
    # 保持后段长度
    out[idx2] <- sign(x[idx2]) * (from + (to - from) / factor + (abs(x[idx2]) - to))
    out
  }
  
  inv <- function(x) {
    out <- x
    cut_pt <- from + (to - from) / factor
    idx1 <- which(!is.na(x) & abs(x) > from & abs(x) < cut_pt)
    idx2 <- which(!is.na(x) & abs(x) >= cut_pt)
    # 还原压缩段
    out[idx1] <- sign(x[idx1]) * (from + (abs(x[idx1]) - from) * factor)
    # 还原后段
    out[idx2] <- sign(x[idx2]) * (to + (abs(x[idx2]) - cut_pt))
    out
  }
  
  trans_new("squash_symmetric", trans, inv)
}

density_umap_triple <- function(
    seu,
    gene1,
    gene2,
    reduction = "umap",
    col_bg   = "#F0F0F0",
    col_g1   = "#155696",  # Plot 1 color (blue)
    col_g2   = "#CCA71B",  # Plot 2 color (yellow)
    col_min  = "#22763F",  # Plot 3 color (green, for example)
    pdf_file = NULL,
    width = 12, height = 4, point_size = 0.2,
    knn_k = 50,
    # New: radial enhancement (center-to-outward fade) switch and strength
    radialize = TRUE,
    radial_scale = 0.90,   # Use the given quantile of positive-point distances as sigma (typ. 0.7–0.95)
    eps = 0
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(grid)
  })
  
  if (is.null(pdf_file) || !nzchar(pdf_file)) {
    pdf_file <- sprintf("DensityUMAP_%s_%s.pdf", gene1, gene2)
  }
  
  stopifnot(reduction %in% Reductions(seu))
  emb <- Embeddings(seu, reduction)
  stopifnot(ncol(emb) >= 2)
  df_base <- data.frame(
    UMAP_1 = emb[,1],
    UMAP_2 = emb[,2],
    cell   = colnames(seu),
    stringsAsFactors = FALSE
  )
  
  rescale01 <- function(x){
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }
  get_expr <- function(g){
    a <- DefaultAssay(seu)
    m <- tryCatch(GetAssayData(seu, assay = a, slot = "data"),
                  error = function(e) GetAssayData(seu, assay = a, slot = "counts"))
    if (!g %in% rownames(m)) stop(sprintf("Gene '%s' is not in the current assay (%s)", g, a))
    as.numeric(m[g, ])
  }
  
  neb_density <- function(features){
    if (!requireNamespace("Nebulosa", quietly = TRUE)) stop("Nebulosa is not installed.")
    ns <- asNamespace("Nebulosa")
    cand <- c("get_density", "calculate_density", ".get_density")
    fn <- NULL
    for (nm in cand) if (exists(nm, where = ns, inherits = FALSE)) { fn <- get(nm, envir = ns); break }
    if (is.null(fn)) stop("Nebulosa internal density function unavailable.")
    tried <- list(
      quote(fn(seu, features = features, reduction = reduction)),
      quote(fn(object = seu, features = features, reduction = reduction)),
      quote(fn(seu, genes = features, reduction = reduction)),
      quote(fn(object = seu, genes = features, reduction = reduction))
    )
    last_err <- NULL
    for (expr in tried) {
      out <- try(eval(expr), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
      last_err <- out
    }
    stop(last_err)
  }
  
  knn_smooth <- function(expr, coords, k = 50){
    if (!requireNamespace("FNN", quietly = TRUE)) {
      stop("Fallback requires the FNN package: install.packages('FNN')")
    }
    nn <- FNN::get.knnx(coords, coords, k)
    sigma <- median(nn$nn.dist[, k], na.rm = TRUE)
    if (!is.finite(sigma) || sigma == 0) sigma <- 1
    w <- exp(-(nn$nn.dist^2) / (2 * sigma^2))
    num <- rowSums(w * expr[nn$nn.index, drop = FALSE])
    den <- rowSums(w)
    val <- num / pmax(den, 1e-8)
    rescale01(val)
  }
  
  # —— Radial enhancement: compute a weighted center and apply isotropic Gaussian decay ——
  radialize_density <- function(d, x, y, scale_q = 0.90, eps = 0){
    d[!is.finite(d)] <- 0
    if (sum(d) <= 0) return(d*0)
    # Weighted centroid
    cx <- sum(x * d) / sum(d)
    cy <- sum(y * d) / sum(d)
    dist <- sqrt((x - cx)^2 + (y - cy)^2)
    # Estimate scale only on positive points (more robust)
    pos <- d > eps
    if (!any(pos)) return(d*0)
    # Sigma = quantile of positive distances (default 0.90; tunable)
    sigma <- as.numeric(quantile(dist[pos], probs = scale_q, na.rm = TRUE))
    if (!is.finite(sigma) || sigma <= 0) sigma <- max(dist[pos], na.rm = TRUE)
    # Radial Gaussian decay
    k <- exp(-(dist^2)/(2*sigma^2))
    rescale01(d * k)
  }
  
  # ------- Compute 0–1 densities d1, d2 for the two genes -------
  d1 <- d2 <- NULL
  dens_mat <- try(neb_density(c(gene1, gene2)), silent = TRUE)
  if (!inherits(dens_mat, "try-error")) {
    d1 <- rescale01(dens_mat[, 1])
    d2 <- rescale01(dens_mat[, 2])
  } else {
    message("[fallback] Using KNN Gaussian smoothing to approximate Nebulosa: ", conditionMessage(attr(dens_mat, "condition")))
    expr1 <- get_expr(gene1)
    expr2 <- get_expr(gene2)
    coords <- as.matrix(df_base[, c("UMAP_1", "UMAP_2")])
    d1 <- knn_smooth(expr1, coords, k = knn_k)
    d2 <- knn_smooth(expr2, coords, k = knn_k)
  }
  
  # —— Optional radialization (more “center-outward” look) ——
  if (radialize) {
    d1 <- radialize_density(d1, df_base$UMAP_1, df_base$UMAP_2, scale_q = radial_scale, eps = eps)
    d2 <- radialize_density(d2, df_base$UMAP_1, df_base$UMAP_2, scale_q = radial_scale, eps = eps)
  }
  
  df <- df_base
  df$d1 <- d1
  df$d2 <- d2
  df$d3 <- pmin(df$d1, df$d2)  # Plot 3 still uses min as intensity
  
  theme_axes_legend <- theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      legend.key.height = unit(18, "pt"),
      legend.key.width  = unit(8,  "pt"),
      legend.background  = element_rect(fill = NA, color = NA),
      legend.box.background = element_rect(fill = NA, color = NA),
      panel.grid = element_blank(),
      axis.line  = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text  = element_text(color = "black", size = 8),
      axis.title = element_text(size = 9),
      plot.background  = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA)
    )
  
  p1 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = d1)) +
    geom_point(size = point_size) +
    scale_color_gradientn(
      colors = c(col_bg, col_g1),
      limits = c(0,1),
      breaks = c(0,1), labels = c("0","1"),
      name   = paste0(gene1, " density")
    ) +
    coord_equal() +
    labs(title = paste0(gene1, " density", if (radialize) " (radialized)" else ""), x = "UMAP_1", y = "UMAP_2") +
    theme_axes_legend
  
  p2 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = d2)) +
    geom_point(size = point_size) +
    scale_color_gradientn(
      colors = c(col_bg, col_g2),
      limits = c(0,1),
      breaks = c(0,1), labels = c("0","1"),
      name   = paste0(gene2, " density")
    ) +
    coord_equal() +
    labs(title = paste0(gene2, " density", if (radialize) " (radialized)" else ""), x = "UMAP_1", y = "UMAP_2") +
    theme_axes_legend
  
  p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = d3)) +
    geom_point(size = point_size) +
    scale_color_gradientn(
      colors = c(col_bg, col_min),
      limits = c(0,1),
      breaks = c(0,1), labels = c("0","1"),
      name   = paste0("overlap (", gene1, ", ", gene2, ")", if (radialize) " (radialized)" else "")
    ) +
    coord_equal() +
    labs(title = paste0("overlap (", gene1, ", ", gene2, ") density", if (radialize) " (radialized)" else ""),
         x = "UMAP_1", y = "UMAP_2") +
    theme_axes_legend
  
  g <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3)
  ggsave(filename = pdf_file, plot = g, width = width, height = height,
         bg = "transparent", useDingbats = FALSE)
  message("Saved: ", normalizePath(pdf_file))
  invisible(g)
}
