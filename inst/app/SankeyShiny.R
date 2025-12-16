# ==========================================
# 1. 基础设置
# ==========================================
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggsankey) # 确保是 davidsjoberg/ggsankey
library(cols4all)
library(colorspace)
library(ggnewscale)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)

# ==========================================
# 2. UI 界面
# ==========================================
ui <- fluidPage(
    # 标题
    titlePanel("Sankey组合图"),
    
    # 侧边布局
    sidebarLayout(
        sidebarPanel(
            width = 3,
            style = "height: 85vh; overflow-y: auto;", 
            
            # --- 1. 数据上传 ---
            h4("1. 数据上传"),
            fileInput("kegg_file", "上传 KEGG 结果 (.csv/.xlsx)"),
            fileInput("diff_file", "上传 差异分析结果 (.csv/.xlsx)"),
            
            # --- 2. 动态列名映射 ---
            h4("2. 列名映射"),
            helpText("系统会自动猜测列名，如有误请修正"),
            uiOutput("col_selectors_kegg"),
            uiOutput("col_selectors_diff"),
            
            # --- 3. 参数设置 ---
            h4("3. 筛选与参数"),
            checkboxInput("convert_id", "ID转换 (Entrez -> Symbol)", value = TRUE),
            uiOutput("pathway_selector"), 
            
            # --- 4. 布局调整 (相对比例) ---
            h4("4. 子图比例调整"),
            sliderInput("width_left", "Bar图 (左) 宽度", 0.1, 2, 1, 0.1),
            sliderInput("width_mid", "Sankey (中) 宽度", 0.1, 2, 0.6, 0.1),
            sliderInput("width_right", "气泡图 (右) 宽度", 0.1, 2, 0.8, 0.1),
            
            # --- 5. 画布尺寸设置 (英寸) ---
            hr(),
            h4("5. 画布尺寸 (英寸/Inch)"),
            helpText("按论文标准设置物理尺寸。预览时按 72DPI 显示，下载 PNG 时按 300DPI 输出。"),
            fluidRow(
                # 【修改】这里单位直接是 inch
                column(6, numericInput("img_width_in", "宽度 (inch)", value = 14, step = 0.5)),
                column(6, numericInput("img_height_in", "高度 (inch)", value = 8, step = 0.5))
            ),
            
            # --- 6. 下载按钮 ---
            div(style="display:flex; gap:10px; flex-wrap:wrap; margin-top: 15px;",
                downloadButton("dl_png", "下载 PNG (300 dpi)"),
                downloadButton("dl_svg", "下载 SVG"),
                downloadButton("dl_pdf", "下载 PDF")
            ),
            
            br(), br()
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel("绘图预览", 
                         # 外层 div 允许滚动，防止大图撑破页面
                         div(style = "overflow: auto; border: 1px solid #eee;",
                             uiOutput("dynamic_plot_ui")
                         )
                ),
                tabPanel("KEGG数据检查", tableOutput("check_kegg")),
                tabPanel("差异数据检查", tableOutput("check_diff"))
            )
        )
    )
)

# ==========================================
# 3. Server 逻辑
# ==========================================
server <- function(input, output, session) {
    
    # --- A. 读取数据 ---
    raw_kegg <- reactive({
        req(input$kegg_file)
        ext <- tools::file_ext(input$kegg_file$name)
        if(ext == "csv") read.csv(input$kegg_file$datapath) else read_excel(input$kegg_file$datapath)
    })
    
    raw_diff <- reactive({
        req(input$diff_file)
        ext <- tools::file_ext(input$diff_file$name)
        if(ext == "csv") read.csv(input$diff_file$datapath) else read_excel(input$diff_file$datapath)
    })
    
    # --- B. 动态列名选择器 ---
    output$col_selectors_kegg <- renderUI({
        req(raw_kegg())
        cols <- colnames(raw_kegg())
        sel_desc <- cols[grep("Desc|Term|Pathway", cols, ignore.case=T)][1]
        sel_gene <- cols[grep("gene|ID", cols, ignore.case=T)][1]
        sel_cnt  <- cols[grep("Count|Rich|Fold", cols, ignore.case=T)][1]
        sel_p    <- cols[grep("pvalue|padj|FDR|qval", cols, ignore.case=T)][1]
        tagList(
            selectInput("c_desc", "Pathway列:", choices=cols, selected=sel_desc),
            selectInput("c_gene", "GeneID列:", choices=cols, selected=sel_gene),
            selectInput("c_cnt",  "Count/X轴列:", choices=cols, selected=sel_cnt), 
            selectInput("c_p",    "Pvalue/颜色列:", choices=cols, selected=sel_p)
        )
    })
    
    output$col_selectors_diff <- renderUI({
        if(is.null(raw_diff())) return(helpText("未上传差异表，将使用随机数据"))
        cols <- colnames(raw_diff())
        sel_gene <- cols[grep("gene|symbol|row", cols, ignore.case=T)][1]
        sel_fc   <- cols[grep("logFC|log2FC|Fold", cols, ignore.case=T)][1]
        sel_p    <- cols[grep("adj|FDR|pvalue", cols, ignore.case=T)][1]
        tagList(
            selectInput("d_gene", "Diff基因列:", choices=cols, selected=sel_gene),
            selectInput("d_fc",   "LogFC列:", choices=cols, selected=sel_fc),
            selectInput("d_p",    "Padj列:", choices=cols, selected=sel_p)
        )
    })
    
    # --- C. 数据清洗 ---
    clean_kegg <- reactive({
        req(raw_kegg(), input$c_desc, input$c_gene, input$c_cnt, input$c_p)
        df <- raw_kegg()
        df$Description <- df[[input$c_desc]]
        df$geneID      <- df[[input$c_gene]]
        df$Count       <- df[[input$c_cnt]] 
        df$pvalue      <- df[[input$c_p]]
        return(df)
    })
    
    # --- D. 通路选择器 ---
    output$pathway_selector <- renderUI({
        req(clean_kegg())
        df <- clean_kegg()
        all_paths <- unique(df$Description)
        selectizeInput("sel_paths", "选择展示的通路:", 
                       choices = all_paths, 
                       selected = head(all_paths, 6), 
                       multiple = TRUE, 
                       options = list(maxItems = 15))
    })
    
    # --- E. 核心绘图逻辑 ---
    plot_obj <- reactive({
        req(clean_kegg(), input$sel_paths)
        
        # 1. 获取数据并筛选
        KEGG_sig <- clean_kegg()
        pathways_select <- input$sel_paths
        KEGG_subset <- KEGG_sig[match(pathways_select, KEGG_sig$Description), , drop = FALSE]
        KEGG_subset <- na.omit(KEGG_subset) 
        KEGG_subset <- KEGG_subset[order(KEGG_subset$Count, decreasing = FALSE), , drop = FALSE]
        KEGG_subset$Pathway <- factor(KEGG_subset$Description, levels = KEGG_subset$Description)
        
        # 2. 构建 sankey 数据
        sankey_gene <- KEGG_subset |>
            as_tibble() |>
            mutate(gene = str_split(geneID, "/")) |>
            tidyr::unnest(cols = gene) |>
            mutate(pathway = paste0(Pathway)) |>
            dplyr::select(gene, pathway, pvalue, Count) 
        
        if(input$convert_id) {
            tryCatch({
                entrez_to_symbol <- bitr(unique(sankey_gene$gene), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
                sankey_gene$gene <- entrez_to_symbol$SYMBOL[match(sankey_gene$gene, entrez_to_symbol$ENTREZID)]
                sankey_gene <- na.omit(sankey_gene)
            }, error = function(e) warning("ID转换失败"))
        }
        
        frame_sankey <- sankey_gene |> make_long(gene, pathway)
        frame_sankey$node <- factor(frame_sankey$node, levels = c(unique(sankey_gene$pathway), unique(sankey_gene$gene)))
        
        # 3. 设置调色
        n_nodes <- length(unique(frame_sankey$node))
        n_flows <- nrow(frame_sankey[!is.na(frame_sankey$next_node), ])
        base_cols <- c4a("kovesi.rainbow_bgyr_35_85_c73", n = n_nodes)
        flow_cols <- c4a("kovesi.rainbow_bgyr_35_85_c73", n = max(1, n_flows))
        my_colors <- lighten(desaturate(base_cols, amount = 0.8), amount = 0.6)
        flow_cols <- lighten(desaturate(flow_cols, amount = 0.8), amount = 0.6)
        flow_cols <- rep(flow_cols, each = 600)
        
        # 4. P1 绘图
        p1 <- ggplot(frame_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
            geom_sankey(flow.fill = flow_cols, flow.color = flow_cols, node.fill = my_colors, width = 0.15, inherit.aes = TRUE) +
            geom_sankey_text(size = 3, color = "black", hjust = 0.5) +
            theme_void() +
            theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
            coord_cartesian(clip = "off") +
            scale_x_discrete(expand = expansion(mult = c(0.1, 0.05)))
        
        # 5. 坐标提取
        pbuilt <- ggplot_build(p1)
        node_data <- tryCatch(pbuilt$data[[2]], error = function(e) NULL)
        if(is.null(node_data)) validate("绘图数据生成失败，请检查是否选择了有效的通路。")
        
        right_nodes <- node_data %>% filter(x == 2) %>% mutate(y = (ymax + ymin) / 2) %>% dplyr::select(label, y, ymin, ymax)
        KEGG_subset$position_y <- right_nodes$y[match(KEGG_subset$Pathway, right_nodes$label)]
        ymin <- min(node_data$ymin)
        ymax <- max(node_data$ymax)
        
        # 6. P2 绘图
        p2 = ggplot(KEGG_subset, aes(x = Count, y = position_y)) +
            scale_color_manual(values = c("KEGG" = "#C5DEBA"), name = "") +
            ggnewscale::new_scale_color() +
            geom_point(aes(y = position_y, color = pvalue, size = Count)) +
            scale_size(range = c(5, 10)) +
            theme_void() +
            scale_color_gradient(low = '#4393C3', high = '#D6604D', name = paste0("-log10(", input$c_p, ")")) +
            labs(x = "", y = "") + 
            theme(
                axis.title = element_text(size = 30, color = "black", family = "sans"),
                axis.text = element_text(size = 20, color = "black", family = "sans"),
                strip.text = element_text(size = 20, face = "bold", color = "black", family = "sans"),
                legend.title = element_text(size = 18, color = "black", family = "sans"),
                legend.text = element_text(size = 15, color = "black", family = "sans"),
                axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                legend.position = "none", axis.ticks.length = unit(0, "pt"), plot.margin = margin(0, 8, 3, 0)
            ) +
            guides(size = "none", alpha = "none") +
            scale_y_continuous(limits = c(ymin, ymax), breaks = KEGG_subset$position_y) +
            scale_x_continuous(limits = c(min(KEGG_subset$Count), max(KEGG_subset$Count) + 2)) +
            geom_segment(aes(x = min(KEGG_subset$Count), xend = min(KEGG_subset$Count), y = min(right_nodes$ymin), yend = max(right_nodes$ymax)), color = "black", size = 0.5) +
            geom_text(data = tibble(y = KEGG_subset$position_y, label = rep("-", length(KEGG_subset$position_y))), aes(x = min(KEGG_subset$Count), y = y, label = label), size = 4, hjust = 1) +
            geom_segment(aes(x = min(KEGG_subset$Count), xend = max(KEGG_subset$Count), y = min(right_nodes$ymin), yend = min(right_nodes$ymin)), color = "black", size = 0.5) +
            geom_text(data = tibble(x = c(min(KEGG_subset$Count), max(KEGG_subset$Count)), y = min(right_nodes$ymin), label = c(round(min(KEGG_subset$Count), 2), round(max(KEGG_subset$Count), 2))), aes(x = x, y = y - 1, label = label), size = 4, vjust = 1, hjust = 0.5) +
            geom_text(data = tibble(x = (min(KEGG_subset$Count) + max(KEGG_subset$Count)) / 2, y = min(right_nodes$ymin), label = input$c_cnt), aes(x = x, y = y, label = label), hjust = 0.5, vjust = 1.2, size = 6, family = "sans")
        
        # 7. P3 绘图
        left_nodes <- node_data %>% filter(x == 1) %>% mutate(y = (ymax + ymin) / 2) %>% dplyr::select(label, y, ymin, ymax)
        gene_list <- as.character(left_nodes$label)
        
        if(is.null(raw_diff())) {
            set.seed(123)
            diffFrame = data.frame(gene = gene_list, log2FC = runif(length(gene_list), 1, 2), padj = runif(length(gene_list), 0.001, 0.05))
            rownames(diffFrame) = diffFrame$gene
        } else {
            df_d <- raw_diff()
            diffFrame <- data.frame(
                gene   = df_d[[input$d_gene]],
                log2FC = df_d[[input$d_fc]],
                padj   = df_d[[input$d_p]]
            )
            diffFrame <- diffFrame[!duplicated(diffFrame$gene), ]
            rownames(diffFrame) <- diffFrame$gene
            diffFrame <- diffFrame[gene_list, ] 
        }
        diffFrame$log2FC[is.na(diffFrame$log2FC)] <- 0
        diffFrame$padj[is.na(diffFrame$padj)] <- 1
        
        left_nodes = cbind(left_nodes, diffFrame)
        
        p3 = ggplot(left_nodes, aes(fill = -log10(padj))) +
            geom_rect(aes(xmin = ymin, xmax = ymax, ymin = 0, ymax = abs(log2FC))) +
            coord_flip() +
            scale_fill_gradient(low = "skyblue", high = "red") +
            scale_y_reverse() +
            theme_classic() +
            theme(
                axis.title = element_blank(), axis.text = element_blank(),
                axis.ticks = element_blank(), axis.line = element_blank(),
                legend.position = "none"
            ) +
            geom_text(aes(x = y, y = 0, label = paste0("log2FC=", round(log2FC, 2), ", padj=", signif(padj, 3))), hjust = 1, vjust = 0.5, size = 3.5)
        
        # 8. 拼接
        layout <- "ACB"
        combined <- p3 + p2 + p1 + 
            patchwork::plot_layout(design = layout, widths = c(input$width_left, input$width_mid, input$width_right))
        
        return(combined)
    })
    
    # --- F. 动态渲染预览图 (逻辑简化) ---
    output$dynamic_plot_ui <- renderUI({
        # 1. 拿英寸
        w_in <- input$img_width_in
        h_in <- input$img_height_in
        
        # 2. 乘72转像素供屏幕显示
        #    使用 ceiling 取整防止小数点导致 CSS 渲染模糊
        w_px <- paste0(ceiling(w_in * 72), "px")
        h_px <- paste0(ceiling(h_in * 72), "px")
        
        # 3. 输出 plotOutput
        plotOutput("main_plot", width = w_px, height = h_px)
    })
    
    output$main_plot <- renderPlot({ plot_obj() })
    
    # --- G. 数据检查 ---
    output$check_kegg <- renderTable({
        req(clean_kegg(), input$sel_paths)
        clean_kegg() %>% filter(Description %in% input$sel_paths) %>% head(10)
    })
    output$check_diff <- renderTable({
        if(is.null(raw_diff())) return(data.frame(Info="未上传差异表，使用随机数据"))
        head(raw_diff(), 10)
    })
    
    # --- H. 智能下载逻辑 (直接使用英寸) ---
    output$dl_png <- downloadHandler(
        filename = "combined_plot.png",
        content = function(file) {
            # 物理尺寸 = 输入的英寸
            # dpi = 300 (高清) -> 像素 = 14 * 300 = 4200px (非常清晰)
            ggsave(file, plot_obj(), width = input$img_width_in, height = input$img_height_in, dpi = 300) 
        }
    )
    
    output$dl_svg <- downloadHandler(
        filename = "combined_plot.svg",
        content = function(file) { 
            ggsave(file, plot_obj(), width = input$img_width_in, height = input$img_height_in) 
        }
    )
    
    output$dl_pdf <- downloadHandler(
        filename = "combined_plot.pdf",
        content = function(file) { 
            ggsave(file, plot_obj(), width = input$img_width_in, height = input$img_height_in, device = "pdf") 
        }
    )
}

shinyApp(ui, server)










