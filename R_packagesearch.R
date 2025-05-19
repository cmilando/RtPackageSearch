library(utils)
library(rvest)
library(httr)

get_packages <- function(term) {
  url <- RSiteSearch(paste0("{",term,"}"),
                            matchesPerPage = 10000)

  Sys.sleep(10)

  # Fetch and parse the HTML
  page <- GET(url, user_agent("Mozilla/5.0")) # avoid bot detection
  html <- read_html(page)

  # Extract all hyperlinks
  links <- html %>%
    html_elements("a") %>%
    html_attr("href")

  # Keep only non-empty and absolute URLs
  links <- links[!is.na(links)]
  links <- unique(links)

  links

  internal_links <- grep("^/CRAN/package", links, value = TRUE)
  x1 <- lapply(internal_links, function(l)
    strsplit(l, "/", fixed = T)[[1]][4])
  y1 <- sort(unique(do.call(c, x1)))

  internal_links <- grep("^/CRAN/refmans", links, value = TRUE)
  x2 <- lapply(internal_links, function(l)
    strsplit(l, "/", fixed = T)[[1]][4])
  y2 <- sort(unique(do.call(c, x2)))


  return(sort(unique(c(y1, y2))))
}

x1 <- get_packages("reproduction number")
x2 <- get_packages("reproduction numbers")
x3 <- get_packages("reproductive number")
x4 <- get_packages("reproductive numbers")
x5 <- get_packages("reproduction rate")
x6 <- get_packages("reproduction rates")

x_all <- sort(unique(c(x1, x2, x3, x4, x5, x6)))
x_all

# ok now, get all package meta data
# package search final

library("pkgsearch")
library("pillar") # nicer data frame printing
cols_to_grab <- c('Package', '')
pkg_data <- lapply(x_all, function(x)
  tail(cran_package_history(x), 1)[, c('Package', 'Title','Date/Publication', 'Description')])

pkg_data[[1]]

pkg_data <- do.call(rbind, pkg_data)
head(pkg_data)
View(pkg_data)


###
install.packages("universe", repos = "https://ropensci.r-universe.dev")
y1 <- universe::global_search(query = '"reproductive number"', limit = 1000L)

y2 <- universe::global_search(query = '"reproductive numbers"', limit = 1000L)

y3 <- universe::global_search(query = '"reproduction number"', limit = 1000L)

y4 <- universe::global_search(query = '"reproduction numbers"', limit = 1000L)

y5 <- universe::global_search(query = '"reproduction rate"', limit = 1000L)

y6 <- universe::global_search(query = '"reproduction rates"', limit = 1000L)

y1$results
# Package, Title, Date/Publication, Description
universe_outputs <- vector("list", 1000)
i_tot = 1
for(i in 1:length(y1$results)) {
  tmp <- data.frame(Package = y1$results[[i]]$Package,
                    Title = y1$results[[i]]$Title,
                    Date = NA,
                    Description = y1$results[[i]]$Description)
  universe_outputs[[i_tot]] = tmp
  i_tot = i_tot + 1
}

for(i in 1:length(y2$results)) {
  tmp <- data.frame(Package = y2$results[[i]]$Package,
                    Title = y2$results[[i]]$Title,
                    Date = NA,
                    Description = y2$results[[i]]$Description)
  universe_outputs[[i_tot]] = tmp
  i_tot = i_tot + 1
}

for(i in 1:length(y3$results)) {
  tmp <- data.frame(Package = y3$results[[i]]$Package,
                    Title = y3$results[[i]]$Title,
                    Date = NA,
                    Description = y3$results[[i]]$Description)
  universe_outputs[[i_tot]] = tmp
  i_tot = i_tot + 1
}

for(i in 1:length(y4$results)) {
  tmp <- data.frame(Package = y4$results[[i]]$Package,
                    Title = y4$results[[i]]$Title,
                    Date = NA,
                    Description = y4$results[[i]]$Description)
  universe_outputs[[i_tot]] = tmp
  i_tot = i_tot + 1
}

for(i in 1:length(y5$results)) {
  tmp <- data.frame(Package = y5$results[[i]]$Package,
                    Title = y5$results[[i]]$Title,
                    Date = NA,
                    Description = y5$results[[i]]$Description)
  universe_outputs[[i_tot]] = tmp
  i_tot = i_tot + 1
}

for(i in 1:length(y6$results)) {
  tmp <- data.frame(Package = y6$results[[i]]$Package,
                    Title = y6$results[[i]]$Title,
                    Date = NA,
                    Description = y6$results[[i]]$Description)
  universe_outputs[[i_tot]] = tmp
  i_tot = i_tot + 1
}

universe_outputs_df <- do.call(rbind, universe_outputs)
universe_outputs_df

r1 <- which(!(universe_outputs_df$Package %in% pkg_data$Package))

universe_outputs_df <- universe_outputs_df[r1, ]
universe_outputs_df
universe_outputs_df <- universe_outputs_df[!duplicated(universe_outputs_df),]
names(pkg_data) = names(universe_outputs_df)
pkg_data2 <- rbind(pkg_data, universe_outputs_df)

library(writexl)
writexl::write_xlsx(pkg_data2, "R_search_v5.xlsx")
