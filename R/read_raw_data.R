library(rio)
library(plyr)
library(tidyverse)

allData <- as.data.frame( import( here::here("data/injection_order_file.xlsx"), sheet = "Organized"))
dataSub <- allData[ , -c( 1:9)] ##Removing the Blanks
rownames( allData) <- allData$compoundId
colnames( allData)
# order the data by injection order (pooled and repeat samples are listed at the end in the original file)
colnames(dataSub)
dataSub.ordered <- dataSub[ , c( 618, 1:25, 619, 26:50, 620, 51:75, 621, 76:100, 622, 101:125, 623, 126:150, 624, 151:175, 625, 176:200, 626, 201:239, 629, 240:264, 630, 265:289, 631, 290:314, 632, 315:339, 633, 340:358, 634, 610:612, 359:380, 635, 381:397, 636, 398:422, 637, 423:447, 638, 448:476, 639, 477:501, 640, 502:526, 641, 527:551, 642, 552:560, 613:617, 561:571, 643, 572:596, 644, 597:609)]
injection_order = data.frame(Run.ID = dataSub.ordered %>% names) %>% mutate(ColPos = 1:length(Run.ID))


raw_data = import(here::here("data/Metab_Metadata.csv"))
raw_batchID = import(here::here("data/fullDatDF.csv")) %>% 
    rename(Run.ID = metabRunNum) %>%   
    select(Run.ID, batch)

raw_data_batch = inner_join(raw_data, raw_batchID, by = "Run.ID")

raw_data %>% dim()
raw_data_batch %>% dim()

raw_data$Run.ID [which(!raw_data$Run.ID %in% raw_batchID$Run.ID)]

raw_data_batch_injection = inner_join(raw_data_batch, injection_order, by = "Run.ID")

raw_data_batch_injection %>% dim()

ldata = raw_data_batch_injection %>%
    select(-V1, -Can.you.read., -Can.you.write., -comments, -Highest.education.level, -datePrepared, -round, -rerun) %>%
pivot_longer(X.10S..Juvenile.hormone.III.acid.diol.isomer:Zearalenone, 
             names_to = "metabolite", values_to = "concentration")

metab = ldata$metabolite[2]

individual_samples = ldata %>%
    filter(! grepl("std_", Run.ID), !grepl("Pool", Run.ID)) %>%
    mutate(Pooled = 0) |>
    mutate(ID = gsub("_2", "", Run.ID)) |>
    arrange(ColPos) |>
    select(ID, ColPos, Pooled, batch, everything())

pooled_samples = ldata |>
    filter(grepl("Pool", Run.ID)) |>
    mutate(ColPos = as.numeric(ColPos), Pooled = 1, ID = Run.ID) |>
    arrange(ColPos) |>
    select(ID, ColPos, Pooled, batch, everything())

column_samples = rbind(individual_samples, pooled_samples) |> arrange(ColPos) |> select(-Run.ID)



ddply(column_samples, .(metabolite),  function(df) {
    missingness = sum(df$concentration < 1e-8) / nrow(df)
    missingness
})

missingness_df =  ldply(unique(column_samples$metabolite), function(metab) {
    df = column_samples %>% filter(metabolite == metab)
    missingness = sum(df$concentration < 1e-8) / nrow(df)
    data.frame(metabolite = metab, missingness = missingness)
})
head(missingness_df)

metab = missingness_df[1,1]
df = column_samples %>% filter(metabolite == metab, batch == 1)
df$concentration

png("test.png", height = 1000, width = 1000)
hist(missingness_df$missingness)
dev.off()

summary(missingness_df)

missingness_df %>% filter(missingness < 0.02) %>% count()

export(missingness_df, here::here("output/missingness.csv"))

export(column_samples, here::here("output/ordered_samples.csv"))
