
library(ggplot2)
library(gridExtra)

time <- rep(1:4, times = 4) 
condition <- rep(c("CMMNE", "MMNE"), each = 8)
measure = rep(c("p1", "p2"), each = 4)


###__________________________________________###
###_________________ PI _____________________###
###__________________________________________###

PI.B = c(0.041, 0.031, 0.021, 0.016,
         0.050, 0.044, 0.041, 0.032,
         0.042, 0.033, 0.022, 0.019,
         0.053, 0.047, 0.040, 0.038)

p1 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = PI.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(pi)) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.8),
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

Mu1.B <- c(2.141, 1.315, 0.819, 0.517, 
           3.131, 1.957, 1.158, 0.946,
           2.488, 1.942, 1.386, 1.521,
           3.418, 2.367, 1.340, 1.077)

p2 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = Mu1.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(mu[1])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.8),
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

Mu2.B <- c(3.510, 3.300, 2.588, 1.939, 
           3.901, 3.621, 2.675, 2.108,
           3.709, 3.731, 3.089, 2.833,
           4.220, 4.039, 3.408, 3.160)

p3 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = Mu2.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(mu[2])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
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

Mu3.B <- c(5.066, 4.540, 4.429, 4.021, 
           5.112, 3.366, 2.819, 2.753,
           5.144, 4.977, 4.736, 4.342,
           5.328, 5.219, 3.693, 3.582)

p4 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = Mu3.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(mu[3])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.81),
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

La1.B <- c(0.644, 0.553, 0.338, 0.264, 
          2.041, 1.736, 1.275, 0.930,
          0.978, 0.711, 0.584, 0.659,
          2.414, 1.894, 1.549, 1.347)

p5 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = La1.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(lambda[1])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.85, 0.8),
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

La2.B <- c(1.600, 1.630, 1.261, 0.943, 
          4.404, 4.359, 3.465, 2.484,
          1.970, 1.906, 1.824, 1.639,
          5.184, 5.076, 4.784, 4.454)

p6 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = La2.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(lambda[2])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
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

La3.B <- c(1.901, 1.533, 1.410, 1.351, 
          3.917, 3.422, 3.515, 2.758,
          2.414, 2.272, 1.790, 1.469,
          4.667, 4.504, 3.857, 3.062)

p7 = ggplot(data.frame(time = time, condition = condition, measure = measure, 
                       score = La3.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(lambda[3])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
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

SI1.B <- c(0.853, 0.563, 0.535, 0.480, 
          0.977, 0.578, 0.518, 0.502,
          2.598, 1.428, 1.162, 0.781,
          2.391, 2.554, 2.410, 3.086,
          3.739, 3.946, 3.428, 4.997,
          7.721, 6.751, 5.050, 7.698)

p8 = ggplot(data.frame(time = rep(1:4, times = 6), 
                       condition = rep(c("CMMNE", "MMNE"), each = 12), 
                       measure = rep(c("p1", "p2", "p3"), each = 4), 
                       score = SI1.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(Sigma[1])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.2, 0.7),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3", "goldenrod2"), name = " ",
                      labels = c(expression(sigma[111]),expression(sigma[112]), 
                                 expression(sigma[122])))

###__________________________________________###
###_________________ SI2 ____________________###
###__________________________________________###

SI2.B <- c(1.668, 1.524, 1.245, 1.172, 
          1.788, 1.983, 1.363, 1.309,
          3.794, 3.214, 2.639, 2.524,
          2.510, 3.135, 4.732, 5.765,
          4.216, 4.676, 6.994, 8.735,
          8.820, 7.907, 11.172, 13.620)

p9 = ggplot(data.frame(time = rep(1:4, times = 6), 
                       condition = rep(c("CMMNE", "MMNE"), each = 12), 
                       measure = rep(c("p1", "p2", "p3"), each = 4), 
                       score = SI2.B), 
            aes(time, score, group = interaction(measure, condition), 
                color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(Sigma[2])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.2, 0.8),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3", "goldenrod2"), name = "",
                      labels = c(expression(sigma[211]),expression(sigma[212]), 
                                 expression(sigma[222])))

###__________________________________________###
###_________________ SI3 ____________________###
###__________________________________________###

SI3.B <- c(1.325, 1.286, 0.939, 0.915, 
          1.710, 1.680, 1.264, 1.190,
          3.854, 3.717, 3.471, 3.429,
          3.835, 4.445, 5.608, 6.883,
          5.330, 5.650, 6.547, 7.463,
          9.809, 10.128, 9.647, 9.816)

p10 = ggplot(data.frame(time = rep(1:4, times = 6), 
                        condition = rep(c("CMMNE", "MMNE"), each = 12), 
                        measure = rep(c("p1", "p2", "p3"), each = 4), 
                        score = SI3.B), 
             aes(time, score, group = interaction(measure, condition), 
                 color = measure, linetype = condition)) +
  geom_line() + geom_point() + ggtitle(expression(Sigma[3])) +
  xlab("n") + ylab("RMSE") + theme_bw() + guides(linetype = FALSE) +
  theme(legend.position = c(0.2, 0.67),
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("100", "200", "500", "1000")) + 
  scale_colour_manual(values = c("brown3", "deepskyblue3", "goldenrod2"), name = "",
                      labels = c(expression(sigma[311]),expression(sigma[312]), 
                                 expression(sigma[322])))


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
                               p9+ rremove("ylab"), p10+ rremove("ylab"), p1 + rremove("ylab"), ncol = 5 , nrow = 2),
                   mylegend, nrow = 2, heights=c(10, 1))
ggsave("fig5MSE.eps", p12, width = 15, height = 8)




