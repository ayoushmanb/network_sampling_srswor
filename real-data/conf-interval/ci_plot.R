# Code to create a side by side plot of CI for comparison
# Load res.jack, res.boot, f in the environment

# pdf("fb_tv_show_cl_ci.pdf")

plot(f-0.005,res.jack[2,], pch = 19, ylim = c(0.2,1.3), xlim = c(0.04,0.32), col = "red", xaxt='n', 
     xlab = "Sampling fraction", ylab = "", cex.axis = 2, cex.lab = 2)
axis(1, at = f, labels = T, cex.axis=2)
points(f+0.005,res.boot[2,], pch = 19, col = "blue")
segments(x0 = f-0.005, y0 = res.jack[1,], x1 = f-0.005, y1= res.jack[3,], col = 'red', lwd = 2)
segments(x0 = f+0.005, y0 = res.boot[1,], x1 = f+0.005, y1= res.boot[3,], col = 'blue', lwd = 2)
abline(h = pop.val, col = "black", lwd = 2)
legend(x = 0.2 , y = 1.3, pch = 19, col = c("red", "blue"), legend = c("Jackknife", "Bootstrap"), cex = 1.5)

# dev.off()