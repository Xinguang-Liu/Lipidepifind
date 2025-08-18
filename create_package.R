if (FALSE) {
  # 激活R文件夹内的所有函数，供测试使用
  devtools::load_all()

  # 加载数据
  library(readxl)
  data1 <- read_xlsx("inst/extdata/Parent lipids.xlsx") %>% as.data.table()
  data2 <- read_xlsx("inst/extdata/Epireactions.xlsx") %>% as.data.table()
  data3 <- read_xlsx("inst/extdata/Allfeatures dataset.xlsx") %>% as.data.table()
  # x_value <- read_excel("inst/extdata/Number of lipid metabolites identified by metabolic reaction (x value).xlsx")
  x_value <- read_excel("inst/extdata/Number of lipid metabolites identified by metabolic reaction-NEW.xlsx")
  result <- read_excel("inst/extdata/result.xlsx")

  saveRDS(data1,"data/data1.rds")
  saveRDS(data2,"data/data2.rds")
  saveRDS(data3,"data/data3.rds")
  saveRDS(x_value,"data/x_value.rds")
  saveRDS(result,"data/result.rds")
  # get_metabolite_mz_predictor(data1,data2,data3,engine = "sqlite")
  lipid_epimetabolite_matching(data1,data2,data3,engine = "sqlite")

  # 为所有函数在man文件夹(如果没有，会创建)下逐一自动建立Rd文档，以及更新NAMESPCAE文档
  devtools::document()

  # 如果没有依赖到别的具有不同版权的第三方包的话，一般选择最为广泛使用的 MIT 即可
  usethis::use_mit_license()


  usethis::use_package("tidyverse", type = "Depends")
  usethis::use_package("data.table", type = "Depends")
  usethis::use_package("dtplyr", type = "Depends")
  usethis::use_package("DBI", type = "Depends")
  usethis::use_package("RSQLite", type = "Depends")
  usethis::use_package("readxl", type = "Depends")

  devtools::check()

  # 存储为rdata格式并使其对用户可用
  usethis::use_data(data1)
  usethis::use_data(data2)
  usethis::use_data(data3)
  usethis::use_data(x_value)
  usethis::use_data(result)

  devtools::document()
  devtools::check()


  data(data1)
}
