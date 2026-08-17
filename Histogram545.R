library(tidyverse)

con <- DBI::dbConnect(
  odbc::odbc(),
  Driver            = "SnowflakeDSIIDriver",   
  Server            = Sys.getenv("SNOWFLAKE_SERVER"),
  warehouse         = Sys.getenv("SNOWFLAKE_WHOUSE"),
  UID               = Sys.getenv("SNOWFLAKE_UID"),
  Authenticator     = "snowflake_jwt",
  PRIV_KEY_FILE     = Sys.getenv("SNOWFLAKE_KEY_PATH"),
  PRIV_KEY_FILE_PWD = Sys.getenv("SNOWFLAKE_KEY_PWD")
)

Batch545EQ <- DBI::dbGetQuery(con, 
                                       "SELECT 
    BATCH_ID
    ,MINOR AS EQ
FROM ROCKIT_DATA_PROD.COMPAC.STG_COMPAC_BATCH
WHERE START_TIME > '2026-01-01 00:00:00.000' 
AND START_TIME <= '2027-01-01 00:00:00.000' 
AND NOT (SIZER_GRADE_NAME IN ('Capture','Captures','Rcy','Capture ','Recycle','Doub','Doubles ','Capt','Leaf','Cap') 
AND NOT SIZE_NAME IN ('US','OS','OOS','OSS'))
AND MINOR BETWEEN 20 AND 110
AND BATCH_ID IN (537,541,545,557)")

DBI::dbDisconnect(con)

#write_csv(Batch545EQ,"BatchEQ.csv")

Batch545EQ |>
  ggplot(aes(EQ)) +
  geom_histogram(binwidth = 0.2) +
  facet_wrap(~BATCH_ID)

#ggsave("HistogramPlot.png", width = 10, height=7)

BatchLevels <- Batch545EQ |> distinct(BATCH_ID) |>
  pull(BATCH_ID)

BatchMoms <- function(Batch,Batch545EQ) {
  
  temp <- Batch545EQ |>
    filter(BATCH_ID == Batch)
  
  outputTemp <- lmom::samlmu(temp$EQ)
  
  as_tibble_row(outputTemp) |> 
    mutate(Batch = Batch, .before = l_1)
  
}

SamLmuAnalysis <- BatchLevels |>
  map(~BatchMoms(.x,Batch545EQ)) |>
  bind_rows()

write_csv(SamLmuAnalysis, "SamLmuAnalysis.csv")  


