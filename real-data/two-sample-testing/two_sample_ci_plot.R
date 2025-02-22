# Code to create a side by side plot of CI for comparison
# Load test_res_jack_mat, test_res_boot_mat, n in the environment

# pdf("test_ci_plot.pdf")

plot(n-25,test_res_jack_mat[2,], pch = 19, ylim = c(0, 20), xlim = c(1800, 3200), col = "red",
     xaxt='n', , yaxt = "n",
     xlab = "n", ylab = "", cex.axis = 2, cex.lab = 2)
axis(1, at = n, labels = T, cex.axis=2)
axis(2, at = c(1,5,10,15,20), labels = T, cex.axis=2)
points(n+25,test_res_boot_mat[2,], pch = 19, col = "blue")
segments(x0 = n-25, y0 = test_res_jack_mat[1,], x1 = n-25, y1= test_res_jack_mat[3,], col = 'red', lwd = 2)
segments(x0 = n+25, y0 = test_res_boot_mat[1,], x1 = n+25, y1= test_res_boot_mat[3,], col = 'blue', lwd = 2)
abline(h = 1, col = "black", lwd = 2)
legend(x = 2600 , y = 19.5, pch = 19, col = c("red", "blue"), legend = c("Jackknife", "Bootstrap"), cex = 1.5)

# dev.off()

