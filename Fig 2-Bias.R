
library(ggplot2)
library(gridExtra)

time <-rep(1:4, times = 4) 
condition <- rep(c("CMMNE", "MMNE"), each = 8)
measure = rep(c("p1", "p2"), each = 4)


###__________________________________________###
###_________________ PI _____________________###
###__________________________________________###

PI.B <- c(-0.012, -0.009, -0.003, -0.002, 
          -0.031, -0.024, -0.020, -0.014,
          -0.014, -0.009, -0.004, -0.006,
          -0.033, -0.026, -0.021,-0.014)

p1 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = PI.B), 
            aes(time, score, group = interaction(measure, condition), 
               color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(pi)) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.2),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = "",
                      labels = c(expression(pi[1]),expression(pi[2])))

###__________________________________________###
###_________________ Mu1 ____________________###
###__________________________________________###

Mu1.B <- c(0.668, 0.316, 0.084, -0.027, 
           -0.965, -0.698, -0.306, -0.134,
          0.869, 0.484, 0.169, 0.233,
          -1.169, -0.662, -0.288,-0.166)

p2 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = Mu1.B), 
            aes(time, score, group = interaction(measure, condition), 
                     color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(mu[1])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.2),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = "",
                      labels = c(expression(mu[11]),expression(mu[12])))

###__________________________________________###
###_________________ Mu2 ____________________###
###__________________________________________###

Mu2.B <- c(0.264, 0.226, 0.364, 0.294, 
           1.050, 0.531, 0.264, 0.096,
           0.299, 0.955, 0.798, 0.851,
           1.399, 0.947, 0.273, 0.205)

p3 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = Mu2.B), 
            aes(time, score, group = interaction(measure, condition), 
                     color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(mu[2])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.8),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = " ",
                      labels = c(expression(mu[21]),expression(mu[22])))

###__________________________________________###
###_________________ Mu3 ____________________###
###__________________________________________###

Mu3.B <- c(3.433, 2.867, 2.659, 2.618, 
           2.174, 1.340, 1.396, 1.303,
           3.468, 3.214, 3.582, 3.481,
           2.224, 1.830, 2.054, 1.697)

p4 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = Mu3.B), 
            aes(time, score, group = interaction(measure, condition), 
                     color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(mu[3])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.42),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = " ",
                      labels = c(expression(mu[31]),expression(mu[32])))


###__________________________________________###
###_________________ La1 ____________________###
###__________________________________________###

La1.B <- c(-0.149, -0.123, -0.056, 0.004, 
           -1.357, -0.998, -0.596, -0.327,
           -0.007, 0.088, 0.191, 0.129,
           -1.042, -0.371, 0.193, 0.101)

p5 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = La1.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(lambda[1])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.2),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = " ",
                      labels = c(expression(lambda[11]),expression(lambda[12])))

###__________________________________________###
###_________________ La2 ____________________###
###__________________________________________###

La2.B <- c(1.201, 1.181, 0.801, 0.430, 
           3.753, 3.543, 2.378, 1.304,
           1.410, 1.291, 1.242, 1.038,
           4.270, 3.955, 3.493, 3.044)

p6 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = La2.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(lambda[2])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.15, 0.5),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = " ",
                      labels = c(expression(lambda[21]),expression(lambda[22])))

###__________________________________________###
###_________________ La3 ____________________###
###__________________________________________###

La3.B <- c(-1.160, -1.031, -0.975, -0.718, 
           -2.788, -2.570, -2.392, -1.499,
           -1.041, -0.900, -0.779, -0.494,
           -2.424, -2.192, -1.819, -1.086)

p7 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = La3.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(lambda[3])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.15, 0.5),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3"), name = " ",
                      labels = c(expression(lambda[31]),expression(lambda[32])))



###__________________________________________###
###_________________ SI1 ____________________###
###__________________________________________###

SI1.B <- c(-0.166, -0.159, -0.150, -0.201, 
           0.346 , 0.210, 0.095, 0.057,
           0.591, 0.272,0.075, -0.067,
           0.373, 0.494, 0.399, 0.696,
           0.816, 0.832, 0.570, 1.131,
           1.957, 1.448, 0.900, 1.668)

p8 = ggplot(data.frame(time = rep(1:4, times = 6), 
                       condition = rep(c("CMMNE", "MMNE"), each = 12), 
                       measure = rep(c("p1", "p2", "p3"), each = 4), 
                       score = SI1.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(Sigma[1])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.6, 0.8),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3", "goldenrod2"), name = " ",
                      labels = c(expression(sigma[1-11]),expression(sigma[1-12]), 
                                 expression(sigma[1-22])))

###__________________________________________###
###_________________ SI2 ____________________###
###__________________________________________###

SI2.B <- c(-0.946, -0.836, -0.544, -0.055, 
           0.840, 0.755, 0.573, 0.545,
           1.877, 1.722, 1.189,  1.081,
           -1.089, 0.399, 1.979, 3.145,
           1.626, 2.079, 3.874, 5.514,
           4.085, 4.052, 6.636, 8.979)

p9 = ggplot(data.frame(time = rep(1:4, times = 6), 
                       condition = rep(c("CMMNE", "MMNE"), each = 12), 
                       measure = rep(c("p1", "p2", "p3"), each = 4), 
                       score = SI2.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(Sigma[2])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.2, 0.8),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3", "goldenrod2"), name = "",
                      labels = c(expression(sigma[2-11]),expression(sigma[2-12]), 
                                 expression(sigma[2-22])))

###__________________________________________###
###_________________ SI3 ____________________###
###__________________________________________###

SI3.B <- c(-0.251, -0.226, -0.110, 0.300, 
           0.627, 0.600, 0.601, 0.562,
           2.151, 2.046, 2.059, 1.826,
           1.461, 1.952, 2.975, 4.231,
           2.402, 2.817, 3.672, 4.721,
           5.647, 6.294, 6.063, 6.382)

p10 = ggplot(data.frame(time = rep(1:4, times = 6), 
                       condition = rep(c("CMMNE", "MMNE"), each = 12), 
                       measure = rep(c("p1", "p2", "p3"), each = 4), 
                       score = SI3.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(Sigma[3])) +
  xlab("n") + ylab("Bias") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.2, 0.67),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3", "goldenrod2"), name = "",
                      labels = c(expression(sigma[3-11]),expression(sigma[3-12]), 
                                 expression(sigma[3-22])))


p11=ggplot(data.frame(time = time, condition = condition, measure = measure, 
                      score = PI.B),  
           mapping = aes(x = 1, y = 1, linetype = condition)) +
  geom_line(size = 1) + theme_void() + theme(legend.position = c(.5, .5), 
                                             legend.text = element_text(size = 13),
                                             legend.direction="horizontal") + theme(
                                               legend.background = element_rect(
                                                 fill = "#FFFFCD", 
                                                 colour = "#FFFFCD", 
                                                 size = 1, linetype="dotted" )) + labs(linetype = "")

library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p11)

library(ggpubr)

p12 = grid.arrange(arrangeGrob(p2, p3+ rremove("ylab"), p4+ rremove("ylab"), 
                         p5+ rremove("ylab"), p6+ rremove("ylab"), p7, p8+ rremove("ylab"), 
                         p9+ rremove("ylab"), p10+ rremove("ylab"), p1 + rremove("ylab"), ncol = 5 , nrow = 2))
ggsave("fig5.eps", p12, width = 15, height = 8)




