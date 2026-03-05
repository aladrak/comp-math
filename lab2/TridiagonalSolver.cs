namespace lab2;

public static class TridiagonalSolver
{
	// Решает систему A x = f, где A — трёхдиагональная матрица.
	// Нижняя диагональ (длина n-1). a[i] — коэффициент при x[i] в (i+1)-м уравнении.
	// Главная диагональ (длина n).
	// Верхняя диагональ (длина n-1). c[i] — коэффициент при x[i+1] в i-м уравнении.
	// Правая часть (длина n).
	// Массив x длины n — решение системы.
	public static double[] Solve(double[] a, double[] b, double[] c, double[] f)
	{
		int n = b.Length;

		// Проверка размеров
		if (a.Length != n - 1)
			throw new ArgumentException($"Длина нижней диагонали должна быть n-1 = {n - 1}, а получена {a.Length}.");
		if (c.Length != n - 1)
			throw new ArgumentException($"Длина верхней диагонали должна быть n-1 = {n - 1}, а получена {c.Length}.");
		if (f.Length != n)
			throw new ArgumentException($"Длина правой части должна быть n = {n}, а получена {f.Length}.");

		// Крайний случай n = 1
		if (n == 1)
		{
			if (Math.Abs(b[0]) < 1e-12)
				throw new DivideByZeroException("Нулевой диагональный элемент в системе размерности 1.");
			return new double[] { f[0] / b[0] };
		}

		// Прогоночные коэффициенты
		double[] alpha = new double[n - 1];
		double[] beta = new double[n - 1];
		double[] x = new double[n];

		// Прямой ход: i = 0
		double denominator = b[0];
		if (Math.Abs(denominator) < 1e-12)
			throw new DivideByZeroException($"Нулевой ведущий элемент на шаге 0.");
		alpha[0] = c[0] / denominator;
		beta[0] = f[0] / denominator;

		// Прямой ход: i = 1 .. n-2
		for (int i = 1; i < n - 1; i++)
		{
			denominator = b[i] - a[i - 1] * alpha[i - 1];
			if (Math.Abs(denominator) < 1e-12)
				throw new DivideByZeroException($"Нулевой ведущий элемент на шаге {i}.");
			alpha[i] = c[i] / denominator;
			beta[i] = (f[i] - a[i - 1] * beta[i - 1]) / denominator;
		}

		// Последняя строка (i = n-1)
		denominator = b[n - 1] - a[n - 2] * alpha[n - 2];
		if (Math.Abs(denominator) < 1e-12)
			throw new DivideByZeroException($"Нулевой ведущий элемент на последнем шаге.");
		x[n - 1] = (f[n - 1] - a[n - 2] * beta[n - 2]) / denominator;

		// Обратный ход
		for (int i = n - 2; i >= 0; i--)
			x[i] = beta[i] - alpha[i] * x[i + 1];

		return x;
	}
}
