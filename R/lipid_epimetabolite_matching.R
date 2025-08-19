#' @title lipid_epimetabolite_matching
#'@description lipid_epimetabolite_matching is a tool to mz matching of lipid epi-metabolites (structure modifications out of classical bio-pathway).
#' @param data1 Parent lipids: We have identified parent lipids from serum,
#' and user can use the results  we provide,
#' or upload your own identified or downloaded parent lipid results from other databases,
#' in the template format of "Parent lipids".
#' @param data2 EpiReactions: The information for the 94 epi-metabolic reactions that may occur in lipids included the reaction type,
#' reaction formula, the change of formula and mass after the reaction,
#' and the matching relationship between the lipid structure and reactions.
#' @param data3 allfeatures dataset: Please upload (differential) ion feature data in the template format of "allfeatures dataset"
#' @param engine "data.table","sqlite"
#'
#' @return result data frame
#' @export
#' @importFrom dplyr %>% mutate across left_join filter sym select
#' @importFrom data.table as.data.table data.table rbindlist
#' @importFrom DBI dbConnect dbWriteTable dbGetQuery dbDisconnect
#' @examples
#' library(readxl)
#' library(Lipidepifind)
#' data_folder <- system.file("extdata", package = "Lipidepifind")
#' data1 <- read_xlsx(file.path(data_folder, "Parent lipids.xlsx")) %>% as.data.table()
#' data2 <- read_xlsx(file.path(data_folder, "Epireactions.xlsx")) %>% as.data.table()
#' data3 <- read_xlsx(file.path(data_folder, "Allfeatures dataset.xlsx")) %>% as.data.table()
#' lipid_epimetabolite_matching(data1,data2,data3,engine = "sqlite")
#'
#' @references 1. Xinguang Liu, Yunfan Zhao, Yang Xie, Yanmin Shi, Jianya Yang, Yan Du, Jinyan Wu, Yitong Gao, Ang Zhang, Jiansheng Li.
#' Characterization of Lipid Epi-metabolites/reaction as Diagnostic Biomarkers for Idiopathic Pulmonary Fibrosis by Lipidepifind, submitted.
lipid_epimetabolite_matching <- function(data1,
                                         data2,
                                         data3,
                                         engine = c("data.table","sqlite")){
  # 加载数据
  data1 <- data1 %>% as.data.table()
  data2 <- data2 %>% as.data.table()
  data3 <- data3 %>% as.data.table()

  # 设置addcution
  addcution1 <- data.table(
    addcution = c(NA, "+NH4", "+H", "-H", "+COOH"),
    value = c(0, 18.0319, 1.0078, -1.0078, 44.9976)
  )


  data2 <- data2 %>%
    # 将数据集2中的加号改为TRUE方便后续筛选
    mutate(across(6:26, ~ . == "+")) %>%
    # 将Mass Difference设置为数字格式
    mutate(`Mass Difference (Da)` = as.numeric(`Mass Difference (Da)`))

  # 设置add
  # 目的是先多对多合并后再筛选符合条件的数据
  data1$add <- 1
  addcution1$add <- 1
  data2$add <- 1
  data3$add <- 1

  # 数据1和addcution多对多合并
  data1 <- data1 %>% left_join(addcution1, by = "add",multiple = "all")

  # 将数据1分为3种，第一种是Double bond小于2考虑compound
  # 第二种Double bond大于等于2
  # 第三种Double bond大于等于3
  data1_1 <- data1 %>% filter(`Double bond FA1` < 2)
  data1_2 <- data1 %>% filter(`Double bond FA1` >= 2)
  data1_3 <- data1 %>% filter(`Double bond FA1` >= 3)



  # 对于第一种需要根据每种compound分别合并数据集2中应该加的行
  class1 <- unique(data1_1$Class)
  a_list <- lapply(seq_along(class1),function(i){
    temp1 <- data1_1 %>% filter(Class == class1[i])
    temp2 <- data2 %>% filter(!!sym(class1[i]) == TRUE)

    temp1 %>%
      left_join(temp2, by = "add",multiple = "all") %>%
      select(Compound, `Molecule weight`, addcution, value, `Mass Difference (Da)`, Reaction) %>%
      mutate(mz = `Molecule weight` + value + `Mass Difference (Da)`)
  })
  # 合并结果生成a_list列表，如用 a_list[[1]] 可以查看compound为FA合并的结果，其中生成的mz就是预测mz
  a <- rbindlist(a_list)

  # 对于第二种需要根据Double bond分别合并数据集2中应该加的行
  temp_b <- data2 %>% filter(`Double bond>=2` == TRUE)
  b <- data1_2 %>%
    left_join(temp_b, by = "add",multiple = "all") %>%
    select(Compound, `Molecule weight`, addcution, value, `Mass Difference (Da)`, Reaction) %>%
    mutate(mz = `Molecule weight` + value + `Mass Difference (Da)`)
  # 合并结果生成b数据集，如用 b 可以查看Double bond>=2合并的结果，其中生成的mz就是预测mz

  # 对于第三种需要根据Double bond分别合并数据集2中应该加的行
  temp_c <- data2 %>% filter(`Double bond>=3` == TRUE)
  c <- data1_3 %>%
    left_join(temp_c, by = "add",multiple = "all") %>%
    select(Compound, `Molecule weight`, addcution, value, `Mass Difference (Da)`, Reaction) %>%
    mutate(mz = `Molecule weight` + value + `Mass Difference (Da)`)

  # 将a, b, c合并为数据集all
  all <- rbind(a,b, c)
  # 此时all为所有数据集1，addcution,数据集2的可能结果

  # 删除无用数据
  rm(temp_b)
  rm(temp_c)
  rm(a)
  rm(b)
  rm(c)

  # 下面用all 和数据集3合并
  # 若直接使用多对多将导致数据量过大，达到2亿多条。
  # 因此在合并前先进行第一步筛选
  # 只有预测mz和数据集3中的m/z相差为1以内的 all数据会被合并
  all$add <- 1

  if (engine == "data.table"){
    # 方法一：lapply data.table ----------------------------------------------------
    all_2_list <- lapply(1:nrow(data3), function(j){
      temp_all1 <- data3[j,]
      temp_value <- temp_all1$`m/z`
      # 此处设置初步筛选条件，只有预测mz和数据集3中的m/z相差为1以内的 all数据会被合并
      temp_all2 <- all %>%
        filter(abs(mz - temp_value) <=1)
      temp_all3 <- merge(temp_all2,temp_all1,by="add")
      print(paste0(j,": ",nrow(temp_all3)))
      temp_all3

    })
    # 合并all_2列表中的所有数据
    all_2 <- rbindlist(all_2_list)
  } else if (engine == "sqlite"){
    # 方法二：使用SQLite数据库条件连接 ---------------------------------------------
    # 创建一个临时SQLite数据库连接
    con <- dbConnect(RSQLite::SQLite(), ":memory:")

    # 将数据表加载到SQLite数据库中
    dbWriteTable(con, "all", all, row.names = FALSE)
    dbWriteTable(con, "data3", data3, row.names = FALSE)

    # 编写SQL合并查询
    query <- '
      SELECT a.*, b.*
      FROM [all] a
      JOIN data3 b ON a.[add] = b.[add]
      WHERE abs(a.mz - b.`m/z`) <= 1
      '

    # 执行SQL查询
    all_2 <- dbGetQuery(con, query)
    all_2 <- all_2[,-length(all_2)] %>% as.data.table()

    # 关闭数据库连接
    dbDisconnect(con)
  }

  # 导出最终数据 -----------------------------------------------------------------

  # 对初步筛选后的数据计算 绝对值ppm，并选择保留ppm小于等于10的数据
  all_3 <- all_2 %>%
    mutate(ppm=abs((mz-`m/z`)/`m/z`)*10^6)%>%
    # 修改下一行设置就可以改变ppm要求，如改为5则变成ppm小于等于5
    filter(ppm<=10)%>%
    select(Compound,`m/z`,Reaction,`Retention time [min]`,Area,addcution,ppm)

  # 导出最终数据
  # write_xlsx(all_3,"result.xlsx")
  return(all_3)
}

