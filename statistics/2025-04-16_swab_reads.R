library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(purrr)
library(ggplot2)
library(rstan)

options(mc.cores = parallel::detectCores())

swab_reads <- read_csv("data/swab_read_counts.csv") %>%
  mutate(date = ymd(paste0("20", date))) %>%
  # assume species assignment is all we care about and that we want deduplicated reads
  summarise(.by = c(date, location, species), viral_reads = sum(dedup_hv), all_reads = first(all_reads)) %>%
  glimpse

pool_data_raw <- read_csv("data/swabs.csv") %>%
  separate_wider_delim(sample, "-", names = c("date", "location"), too_many = "drop") %>%
  mutate(date = ymd(paste0("20", date)))

pool_positivity <- pool_data_raw %>%
  pivot_longer(-c(date, location, pool_size), values_to = "positivity") %>%
  # only get the species-level positivity
  semi_join(swab_reads, by = join_by(name == species)) %>%
  rename(species = name)

# Check positivity == viral_reads < 0
pool_positivity %>%
  left_join(swab_reads, by = join_by(date, location, species)) %>%
  mutate(expected = (positivity == 0 & is.na(viral_reads)) | (positivity == 1 & !is.nan(viral_reads))) %>%
  pull(expected) %>%
  all

pool_data_raw %>%
  select(date, location, pool_size) %>%
  left_join(swab_reads %>% select(date, location, all_reads) %>% distinct)
  

pool_data <- swab_reads %>%
  select(date, location, all_reads) %>%
  distinct %>%
  left_join(select(pool_data_raw, date, location, pool_size)) %>%
  mutate(reads_per_swab = all_reads / pool_size) 

pool_data %>%
  print(n = Inf)
# The reads per sample vary over 2 oom! min = 83, max = 20,000

swab_reads %>%
  mutate(single_digits = viral_reads <= 10) %>%
  summarise(single_digit_frac = mean(single_digits), .by = species)

df <- swab_reads %>%
  left_join(pool_data) %>%
  mutate(ra_per_swab = viral_reads / reads_per_swab) %>%
  arrange(species, date) %>%
  print(n = Inf)

df %>%
  ggplot(mapping = aes(ra_per_swab, species)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 0.02))

df %>%
  ggplot(mapping = aes(date, ra_per_swab)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 0.02))

df %>%
  ggplot(mapping = aes(all_reads, viral_reads)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 200))

df %>%
  ggplot(mapping = aes(reads_per_swab, viral_reads)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 0.002) +
  scale_y_continuous(limits = c(0, 200))

df %>%
  ggplot(mapping = aes(ra_per_swab, viral_reads)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 0.02)) +
  scale_y_continuous(limits = c(0, 200))
