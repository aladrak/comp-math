# x_data <- c(0, 1, 2, 3, 4, 5)
# y_data <- c(1.2, 2.8, 4.1, 5.0, 4.9, 3.1)

x_data <- c(5, 7, 9, 11, 12, -3, 10)
y_data <- c(23, 47, 79, 119, 142, 7, 98)

# Вектор точек для гладкого построения графиков
x_plot <- seq(min(x_data), max(x_data), length.out = 500)

# 2. ИНТЕРПОЛЯЦИЯ ПО ФОРМУЛЕ ЛАГРАНЖА
lagrange_interp <- function(x_nodes, y_nodes, x) {
  n <- length(x_nodes)
  y <- numeric(length(x))
  for (i in seq_along(x)) {
    s <- 0
    for (j in 1:n) {
      l <- 1
      for (k in 1:n) {
        if (k != j) {
          l <- l * (x[i] - x_nodes[k]) / (x_nodes[j] - x_nodes[k])
        }
      }
      s <- s + y_nodes[j] * l
    }
    y[i] <- s
  }
  return(y)
}

y_lagrange <- lagrange_interp(x_data, y_data, x_plot)

# 3. ИНТЕРПОЛЯЦИЯ ПО ФОРМУЛЕ НЬЮТОНА
newton_interp <- function(x_nodes, y_nodes, x) {
  n <- length(x_nodes)
  dd <- matrix(0, n, n)
  dd[, 1] <- y_nodes
  for (j in 2:n) {
    for (i in 1:(n - j + 1)) {
      dd[i, j] <- (dd[i + 1, j - 1] - dd[i, j - 1]) / (x_nodes[i + j - 1] - x_nodes[i])
    }
  }
  a <- dd[1, ]
  
  y <- numeric(length(x))
  for (i in seq_along(x)) {
    res <- a[1]
    prod_term <- 1
    for (j in 2:n) {
      prod_term <- prod_term * (x[i] - x_nodes[j - 1])
      res <- res + a[j] * prod_term
    }
    y[i] <- res
  }
  return(y)
}

y_newton <- newton_interp(x_data, y_data, x_plot)

# 4. СГЛАЖИВАЮЩИЕ МНОГОЧЛЕНЫ (МНК, степени 1-3)
fit1 <- lm(y_data ~ x_data)
fit2 <- lm(y_data ~ x_data + I(x_data^2))
fit3 <- lm(y_data ~ x_data + I(x_data^2) + I(x_data^3))

nd <- data.frame(x_data = x_plot)
y_mnk1 <- predict(fit1, newdata = nd)
y_mnk2 <- predict(fit2, newdata = nd)
y_mnk3 <- predict(fit3, newdata = nd)

# 5. ПРОИЗВОЛЬНЫЙ МНОГОЧЛЕН 4-Й СТЕПЕНИ
# Коэффициенты задаются от свободного члена до старшей степени: 
# P(x) = c[1] + c[2]*x + c[3]*x^2 + c[4]*x^3 + c[5]*x^4
coefs4 <- c(3, 2, 0, 0, 0)

poly4_eval <- function(coefs, x) {
  sapply(x, function(xi) sum(coefs * xi^(0:4)))
}

y_poly4 <- poly4_eval(coefs4, x_plot)

# 6. ПОСТРОЕНИЕ ГРАФИКОВ
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 1, 0))

# 6.1 Лагранж
plot(x_plot, y_lagrange, type = "l", col = "steelblue", lwd = 2,
     main = "Интерполяция Лагранжа", xlab = "x", ylab = "y")
points(x_data, y_data, col = "red", pch = 19, cex = 1.2)
grid()

# 6.2 Ньютон
plot(x_plot, y_newton, type = "l", col = "darkgreen", lwd = 2,
     main = "Интерполяция Ньютона", xlab = "x", ylab = "y")
points(x_data, y_data, col = "red", pch = 19, cex = 1.2)
grid()

# 6.3 МНК (1, 2, 3 степени)
plot(x_plot, y_mnk1, type = "l", col = "blue", lwd = 2, 
     ylim = range(c(y_data, y_mnk1, y_mnk2, y_mnk3)),
     main = "Сглаживание МНК", xlab = "x", ylab = "y")
lines(x_plot, y_mnk2, col = "orange", lwd = 2, lty = 2)
lines(x_plot, y_mnk3, col = "purple", lwd = 2, lty = 3)
points(x_data, y_data, col = "red", pch = 19, cex = 1.2)
legend("topright", legend = c("1 степень", "2 степень", "3 степень", "Узлы"),
       col = c("blue", "orange", "purple", "red"), lty = c(1, 2, 3, NA),
       pch = c(NA, NA, NA, 19), lwd = 2, bg = "white")
grid()

# 6.4 Многочлен 4-й степени
plot(x_plot, y_poly4, type = "l", col = "black", lwd = 2,
     main = "Многочлен 4-й степени (по коэффициентам)", xlab = "x", ylab = "y")
grid()

par(mfrow = c(1, 1))