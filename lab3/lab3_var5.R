# x <- c(5, 7, 9, 11, 12, -3, 10)
# y <- c(3, -2, -2, 4, 15, 4, 1)

x <- c(5, 7, 9, 11, 12, -3, 10)
y <- c(23, 47, 79, 119, 142, 7, 98)

# 1. РАСЧЕТ КОЭФФИЦИЕНТОВ СГЛАЖИВАЮЩИХ МНОГОЧЛЕНОВ 
#    (МЕТОД НАИМЕНЬШИХ КВАДРАТОВ)
# Многочлен степени 1
poly1 <- lm(y ~ poly(x, 1, raw = TRUE))
coef1 <- coef(poly1)
cat("Коэффициенты многочлена степени 1 (МНК):\n")
cat("a0 =", coef1[1], "\n")
cat("a1 =", coef1[2], "\n\n")

# Многочлен степени 2
poly2 <- lm(y ~ poly(x, 2, raw = TRUE))
coef2 <- coef(poly2)
cat("Коэффициенты многочлена степени 2 (МНК):\n")
cat("a0 =", coef2[1], "\n")
cat("a1 =", coef2[2], "\n")
cat("a2 =", coef2[3], "\n\n")

# Многочлен степени 3
poly3 <- lm(y ~ poly(x, 3, raw = TRUE))
coef3 <- coef(poly3)
cat("Коэффициенты многочлена степени 3 (МНК):\n")
cat("a0 =", coef3[1], "\n")
cat("a1 =", coef3[2], "\n")
cat("a2 =", coef3[3], "\n")
cat("a3 =", coef3[4], "\n\n")

# 2. РАСЧЕТ КОЭФФИЦИЕНТОВ ИНТЕРПОЛЯЦИОННОГО МНОГОЧЛЕНА
#    (РЕШЕНИЕ СИСТЕМЫ УРАВНЕНИЙ 3.1 - МАТРИЦА ВАНДЕРМОНДА)
n <- length(x) - 1

# Создание матрицы Вандермонда
V <- matrix(0, nrow = length(x), ncol = length(x))
for (i in seq_along(x)) {
  for (j in seq_along(x)) {
    V[i, j] <- x[i]^(j-1)
  }
}
cat("Матрица Вандермонда:\n")
print(V)
cat("\n")

# Решение V * a = y
coef_interp <- solve(V, y)

cat("Коэффициенты интерполяционного многочлена (степень 4):\n")
for (i in seq_along(coef_interp)) {
  cat("a", i-1, " = ", coef_interp[i], "\n", sep = "")
}
cat("\n")

# 3. ПОСТРОЕНИЕ ГРАФИКОВ
# Создаем последовательность точек для плавных кривых
x_seq <- seq(min(x) - 1, max(x) + 1, length.out = 200)

# Функция для вычисления значения многочлена
eval_poly <- function(x_vals, coefs) {
  result <- rep(0, length(x_vals))
  for (i in 1:length(coefs)) {
    result <- result + coefs[i] * x_vals^(i-1)
  }
  return(result)
}

y_poly1 <- eval_poly(x_seq, coef1)
y_poly2 <- eval_poly(x_seq, coef2)
y_poly3 <- eval_poly(x_seq, coef3)
y_interp <- eval_poly(x_seq, coef_interp)

png("./lab3/polynomials_comparison.png", width = 1200, height = 800)
par(mfrow = c(2, 2))

# График 1: Линейная аппроксимация (степень 1)
plot(x, y, pch = 19, col = "red", cex = 2, 
     main = "Линейная аппроксимация (степень 1)",
     xlab = "x", ylab = "y", xlim = range(x_seq), ylim = range(c(y, y_poly1)))
grid()
lines(x_seq, y_poly1, col = "blue", lwd = 2)
points(x, y, pch = 19, col = "red", cex = 2)
legend("topleft", legend = c("Исходные данные", "Многочлен степени 1"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))

# График 2: Квадратичная аппроксимация (степень 2)
plot(x, y, pch = 19, col = "red", cex = 2,
     main = "Квадратичная аппроксимация (степень 2)",
     xlab = "x", ylab = "y", xlim = range(x_seq), ylim = range(c(y, y_poly2)))
grid()
lines(x_seq, y_poly2, col = "blue", lwd = 2)
points(x, y, pch = 19, col = "red", cex = 2)
legend("topleft", legend = c("Исходные данные", "Многочлен степени 2"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))

# График 3: Кубическая аппроксимация (степень 3)
plot(x, y, pch = 19, col = "red", cex = 2,
     main = "Кубическая аппроксимация (степень 3)",
     xlab = "x", ylab = "y", xlim = range(x_seq), ylim = range(c(y, y_poly3)))
grid()
lines(x_seq, y_poly3, col = "blue", lwd = 2)
points(x, y, pch = 19, col = "red", cex = 2)
legend("topleft", legend = c("Исходные данные", "Многочлен степени 3"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))

# График 4: Интерполяционный многочлен (степень 4)
plot(x, y, pch = 19, col = "red", cex = 2,
     main = "Интерполяционный многочлен (степень 4)",
     xlab = "x", ylab = "y", xlim = range(x_seq), ylim = range(c(y, y_interp)))
grid()
lines(x_seq, y_interp, col = "blue", lwd = 2)
points(x, y, pch = 19, col = "red", cex = 2)
legend("topleft", legend = c("Исходные данные", "Интерполяционный многочлен"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
dev.off()

# ГРАФИК 2: Сравнение всех многочленов
png("./lab3/all_polynomials.png", width = 1000, height = 700)
plot(x, y, pch = 19, col = "black", cex = 2,
     main = "Сравнение всех многочленов",
     xlab = "x", ylab = "y", xlim = range(x_seq), 
     ylim = range(c(y, y_poly1, y_poly2, y_poly3
    #  , y_interp
     )))
grid()
lines(x_seq, y_poly1, col = "blue", lwd = 2, lty = 1)
lines(x_seq, y_poly2, col = "green", lwd = 2, lty = 2)
lines(x_seq, y_poly3, col = "orange", lwd = 2, lty = 3)
# lines(x_seq, y_interp, col = "red", lwd = 2, lty = 4)
points(x, y, pch = 19, col = "black", cex = 2)
legend("topleft",
       legend = c("Исходные данные",
                  "МНК: степень 1",
                  "МНК: степень 2",
                  "МНК: степень 3",
                  "Интерполяция: степень 4"),
       col = c("black", "blue", "green", "orange", "red"),
       pch = c(19, NA, NA, NA, NA),
       lty = c(NA, 1, 2, 3, 4),
       lwd = c(NA, 2, 2, 2, 2))
dev.off()

cat("ИТОГОВЫЕ РЕЗУЛЬТАТЫ:\n\n")
cat("1. МНОГОЧЛЕН СТЕПЕНИ 1 (МНК):\n")
cat(sprintf("   P1(x) = %.4f + %.4f*x\n\n", coef1[1], coef1[2]))
cat("2. МНОГОЧЛЕН СТЕПЕНИ 2 (МНК):\n")
cat(sprintf("   P2(x) = %.4f + %.4f*x + %.4f*x^2\n\n", coef2[1], coef2[2], coef2[3]))
cat("3. МНОГОЧЛЕН СТЕПЕНИ 3 (МНК):\n")
cat(sprintf("   P3(x) = %.4f + %.4f*x + %.4f*x^2 + %.4f*x^3\n\n", 
            coef3[1], coef3[2], coef3[3], coef3[4]))
cat("4. ИНТЕРПОЛЯЦИОННЫЙ МНОГОЧЛЕН (степень 4):\n")
cat(sprintf("   P4(x) = %.4f", coef_interp[1]))
for (i in 2:length(coef_interp)) {
  cat(sprintf(" + %.4f*x^%d", coef_interp[i], i-1))
}

cat("\n\nПРОВЕРКА ИНТЕРПОЛЯЦИИ:\n")
y_check <- eval_poly(x, coef_interp)
for (i in seq_along(x)) {
  cat(sprintf("   P4(%g) = %.4f (исходное: %g)\n", x[i], y_check[i], y[i]))
}