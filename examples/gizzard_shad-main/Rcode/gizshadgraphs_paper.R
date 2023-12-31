### Graph And Analysis of Gizzard Shad Model
source("gizshadmodel.R")

#######################
## Graphing Fun Time!
#######################

## Setup ##
N <- 400 # number of size classes
l_shad <- 0.00   # lower size limit in mm
u_shad <- 500.0    # upper size limit in mm - we want this to be
                   # larger than L-infty
delta_z <- (u_shad - l_shad) / N
zmesh <-  l_shad + ((1:N) - 1 / 2) * (u_shad - l_shad) / N
tf <-100 # number of years

# Initial length distribution
n <- matrix(0, length(zmesh), tf)
n0_total <-1487.929
# n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5*m_par$grow_max, sd = 30)
    #normal like LTRM 1994
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Note: sum(n[,1])*delta_z = n0_total

## Egg Function vs. Data
plot_df <- data_frame(z = zmesh, eggs =  eggs_z(zmesh, m_par)/1000 )
ggplot(data = plot_df,
       aes( x = z, y = eggs)) +
  geom_line(color = "blue", size = 1)+
  labs(x = "length (in mm)",
       y = "eggs (in thousands)",
       title = "Eggs Produced") +
  #     subtitle = "Gizzard Shad")  
  scale_x_continuous(limits = c(0,u_shad), breaks = seq(0,500,100))+
  scale_y_continuous(limits = c(0,700), breaks = seq(0,700,100))+
  geom_point(data = egg_size_data,
             aes(x = x, y = EggLengthData),
             color = "black")+
  theme_classic()+
  theme(text = element_text(size=20),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/eggs.png")

## Survival Probability of Age-0 vs data
dmesh <- seq( from = min(Michaletz_data$density), 
              to =max(Michaletz_data$density), 
              length.out = length(Michaletz_data$density) )

plot_df <- data_frame(x = dmesh, prob = surv_density(dmesh, m_par) )
ggplot(data = plot_df,
       aes( x = x, y = prob))+
  geom_line(color = "blue", size = 1)+
  labs(x = "density (age-0 per 1000 m^3)",
       y = "probability of survival",
       title = "Survival Probability for Age-0")+
  scale_x_continuous(limits = c(0,max(Michaletz_data$density)), 
                     breaks = seq(0,max(Michaletz_data$density),200))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.1))+
  geom_point(data = Michaletz_data, 
             aes(x = density, y = survival/100))+
  theme_classic()+
  theme(text = element_text(size=20),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/age0surv.png")

# Dynamical System
for (i in 1:(tf - 1)) {
  k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                m_par)) * delta_z
  n[, i + 1] <- k_iter %*% n[, i]
}

# Population size vs time
n_total <- rep(0, tf)
for (i in 1:tf) {
  n_total[i] <- sum(n[, i]) * delta_z
}

plot_df <- data.frame(year = 1:tf, total = n_total)
ggplot(plot_df,
       aes(x = year, y = total)) + 
  geom_line() +
  labs(x = "years",
       y = "total") +
  labs(x = "time (in years)",
       y = "total density",
       color = "Legend") +
#  theme_bw() + 
   theme_classic() +  
  theme(legend.position = c(0.8, 0.4))+
  theme(text = element_text(size=20),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/ntotal.png")

### Compute period of n_total
# Make tibble of many periods at/near stable state
many_period <- tibble(
  t = (tf-30):tf,
  n_tot = n_total[(tf-30):tf]
)
# fit sine curve to many_period
n_total_period <- nls(
  n_tot~(mean(many_period$n_tot)+(max(many_period$n_tot)-mean(many_period$n_tot))*sin((2*pi/per)*(t-shift))),
  data = many_period,
  start = list(
    per = 9,
    shift = 5))

# show_years <- tf -11 + seq(from = 1, 
#                        length.out = 5,
#                        by  = 2)
# n_freq <- sweep(n, 2, colSums(n),  FUN = "/")
# plot_df <- data.frame(z = zmesh, year = n_freq[,show_years])
# plot_df <- melt(plot_df, id.vars = 'z', variable.name = 'year')
# ggplot(plot_df,
#        aes(z, value)) + 
#   geom_line(aes(color = year)) + 
#   labs(x = "length (in mm)",
#        y = "relative frequency",
#        title = "n(z,t)/total",
#        color = "Legend") +
#   scale_color_manual(
#     values = rainbow(5),
#     labels = c("Year 1", "Year 3",
#                "Year 5","Year 7", "Year 9"))+
#   theme_bw() +  
#   theme(legend.position = c(0.85, 0.5))+
#   theme(text = element_text(size=16),
#         aspect.ratio = .7)
# ggsave("~/Documents/research/gizzard_shad/figures/period.png")
# 
# show_years <- tf -10 + seq(from = 1, 
#                           length.out = 5,
#                           by  = 1)
# n_freq <- sweep(n, 2, colSums(n),  FUN = "/")
# plot_df <- data.frame(z = zmesh, year = n_freq[,show_years])
# plot_df <- melt(plot_df, id.vars = 'z', variable.name = 'year')
# ggplot(plot_df,
#        aes(z, value)) + 
#   geom_line(aes(color = year)) + 
#   labs(x = "length (in mm)",
#        y = "relative frequency",
#        title = "n(z,t)/total",
#        color = "Legend") +
#   scale_color_manual(
#     values = rainbow(5),
#     labels = c("Year 1", "Year 2",
#                "Year 3","Year 4", "Year 5"))+
#   theme_bw() +  
#   theme(legend.position = c(0.85, 0.5))+
#   theme(text = element_text(size=16),
#         aspect.ratio = .7)
# ggsave("~/Documents/research/gizzard_shad/figures/period_small.png")
# show_years <- tf -6 + seq(from = 1, 
#                           length.out = 5,
#                           by  = 1)
# n_freq <- sweep(n, 2, colSums(n),  FUN = "/")
# plot_df <- data.frame(z = zmesh, year = n_freq[,show_years])
# plot_df <- melt(plot_df, id.vars = 'z', variable.name = 'year')
# ggplot(plot_df,
#        aes(z, value)) + 
#   geom_line(aes(color = year)) + 
#   labs(x = "length (in mm)",
#        y = "relative frequency",
#        title = "n(z,t)/total",
#        color = "Legend") +
#   scale_color_manual(
#     values = rainbow(5),
#     labels = c("Year 5", "Year 6",
#                "Year 7","Year 8", "Year 9"))+
#   theme_bw() +  
#   theme(legend.position = c(0.85, 0.5))+
#   theme(text = element_text(size=16),
#         aspect.ratio = .7)
# ggsave("~/Documents/research/gizzard_shad/figures/period_large.png")
# 

show_years <- tf -13 + seq(from = 1, 
                           length.out = 9,
                           by  = 1)
n_freq <- sweep(n, 2, colSums(n),  FUN = "/")
plot_df <- data.frame(z = zmesh, Year = n_freq[,show_years])
plot_df <- melt(plot_df, id.vars = 'z', variable.name = 'Year')
ggplot(plot_df,
       aes(z, value)) + 
  geom_line() + 
  labs(x = "length (in mm)",
       y = "relative frequency",
       color = "Legend") +
 # scale_color_manual(
#    values = rainbow(5),
#    labels = c("Year 1", "Year 2",
#               "Year 3","Year 4", "Year 5"))+
  facet_wrap(~ Year, nrow = 3)+
  scale_x_continuous(limits = c(0, u_shad), breaks = seq(0, 500, 200),
                     expand = c(0, 0)) +
#  scale_y_continuous(limits = c(0, 0.0), breaks = seq(0, 0.015, 0.005),
#                     expand = c(0, 0))+
  theme_bw() +  
#  theme_classic()+
  theme(legend.position = c(0.85, 0.5))+
  theme(text = element_text(size=20),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/period_facet.png")


### Age-0 survival vs time
surv_t <- rep(0, times = tf)
for (i in 1:tf) {
  surv_t[i] <- surv_age0(n = n[, i], z = zmesh, m_par)
}

plot_df <- tibble(time_years = 1:100, prob = surv_t[1:100])
ggplot(data = plot_df,
       aes(x = time_years, y = prob)) +
  geom_line() +
  labs(x = "time (in years)",
       y = "probability of survival")+
  scale_x_continuous(limits = c(5, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, .04), breaks = seq(0, 0.2, 0.005)) +
  #  theme_bw() +  
  theme_classic()+
  theme(text = element_text(size = 20),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/Figure2a.pdf")

#############################
########
ltrm_gzsd <- read_csv("ltrm_fish_data.csv")
# Remove length 0 and NA
ltrm_gzsd <- ltrm_gzsd[!is.na(ltrm_gzsd$length) & (ltrm_gzsd$length > 0) &
                         !is.na(ltrm_gzsd$fdate), ]
# Convert date into new format
# Then pull year and add it as a new column
ltrm_gzsd$fdate <- as.Date(ltrm_gzsd$fdate, "%m/%d/%Y")
ltrm_gzsd <- ltrm_gzsd %>% mutate(year = year(fdate))
ltrm_gzsd <- ltrm_gzsd %>%
  filter(year != 2107)

# Facet graph of LTRM length frequency in main channel
ltrm_gzsd %>%
  filter(year < 2107) %>%
  filter(pool %in% c("04", "08", "13", "26", "OR")) %>%
  ggplot(., aes(x = length)) +
  geom_histogram(aes(y = ..density..), bins = 30) +
  facet_wrap(~ year, nrow = 4) +
  labs(x = "length (in mm)",
       y = "length frequency",
       subtitle = "Upper Mississippi River - Main Channel, 1989-2020") +
  scale_x_continuous(limits = c(0, u_shad), breaks = seq(0, 500, 200),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.015), breaks = seq(0, 0.015, 0.005),
                     expand = c(0, 0))+
    theme_bw() +  
  #theme_classic() +
  theme(text = element_text(size=14),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/LTRMmain.png")


# Facet graph of LTRM length frequency in main channel
ltrm_gzsd %>%
  filter(year < 2107) %>%
  filter(pool %in% c("LG")) %>%
  ggplot(., aes(x = length)) +
  geom_histogram(aes(y = ..density..), bins = 30) +
  facet_wrap(~ year, ncol = 9) +
  labs(x = "length (in mm)",
       y = "length frequency",
       subtitle = "Upper Mississippi River - La Grange, 1992-2020") +
  scale_x_continuous(limits = c(0, u_shad), breaks = seq(0, 500, 200),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.014), breaks = seq(0, 0.015, 0.005),
                     expand = c(0, 0))+ 
  theme_bw() +
  theme(text = element_text(size=12),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/LTRMlg.png")


## Average of periodic orbit compared with LTRMP data from La Grange
#### Now Check with La Grange
n <- matrix(0, length(zmesh), tf)
n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5 * m_par$grow_max, sd = 30)
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Dynamical System
for (i in 1:(tf - 1)) {
  k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                m_par)) * delta_z
  n[, i + 1] <- k_iter %*% n[, i]
}

year_start <- tf-8
year_end <- tf
plot_average <- tibble(z = zmesh, 
                       year = rep("mean",N),
                       n_freq = (1/9)*rowSums(n_freq[,year_start:year_end]) )
ltrm_gzsd %>%
  filter(year < 2107) %>%
  filter(pool %in% c("LG")) %>%
ggplot(., aes(x = length)) +
  geom_histogram(aes(y = ..density..), bins = 50)+
  #  geom_density(aes(x = length)) +
  geom_line(data = plot_average, aes(x = z, y = n_freq/delta_z))+
  labs(x = "length (in mm)",
       y = "relative frequency",
       color = "Legend") +
  scale_color_manual(breaks=c("a","b"))+
  #theme_bw()+
  theme_classic()+
  theme(text = element_text(size=20),
        aspect.ratio = .7)
ggsave("~/Documents/research/gizzard_shad/figures/lagrange.pdf")

