##### Simulate 1 million points, write into "cvxReg_input.csv"
##### Run convex regression C code, which reads "cvxReg_input.csv" and write fitted curve into "cvxReg_output_1.csv"
##### Plot the fitted curve

### Simulation
set.seed(2017)

n = 1e6

x = seq(0, 1, length = n)
y = exp(-3*x) + 0.1*rnorm(n, 0, 1)
w = rep(1, n)
DataFrame = data.frame(y, x, w)
write.csv(DataFrame, file = "cvxReg_input.csv", quote = FALSE, row.names = FALSE)

### Run C code ("CvxReg.c")

### Plot fitted curve
output_1 = read.csv("cvxReg_output_1.csv")

pdf('CvxSimuRes.pdf', height = 10, width =  15)
par(mfrow = c(1,2), mar = c(10, 10, 3, 1))
plot(x, y, ylim = range(y), ylab = "y", cex.lab = 2)
lines(x, exp(-3*x), col = "red", lwd = 2)
lines(x, output_1[,1], col = "blue", lwd = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
legend("bottomright", col = c("red", "blue"), lwd = 2, lty = c(1,1), legend = c("Original", "Fitted"))

plot(x, output_1[,1], ylim = range(y), ylab = "fitted y", col = "blue", lwd = 2, type = 'l', cex.lab = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
dev.off()

