namespace lab2;
public enum NormType
{
	L1,      // сумма модулей / максимум суммы по столбцам
	L2,      // евклидова / Фробениуса (корень из суммы квадратов)
	LInf     // максимум модуля / максимум суммы по строкам
}
public static class SimpleIteration
{
	/// Решает систему Ax = b методом простых итераций.
	public static double[] Solve(double[,] A, double[] b, double epsilon, int maxIterations = 1000, NormType normType = NormType.LInf)
	{
		int n = A.GetLength(0);
		if (n != A.GetLength(1))
			throw new ArgumentException("Матрица A должна быть квадратной.");
		if (b.Length != n)
			throw new ArgumentException("Размер вектора b не совпадает с размерностью матрицы.");

		// Построение матрицы α и вектора β
		var (alpha, beta) = BuildAlphaBeta(A, b);

		// Проверка достаточного условия сходимости (норма α < 1)
		double normAlpha = Norm(alpha, normType);
		if (normAlpha >= 1.0)
		{
			Console.WriteLine($"Предупреждение: ||α|| = {normAlpha:F4} ≥ 1. Сходимость не гарантирована.");
		}

		// Начальное приближение (β)
		double[] xPrev = (double[])beta.Clone();
		double[] xCurr = new double[n];
		int iteration = 0;

		// Коэффициент для апостериорной оценки погрешности
		double factor = normAlpha / (1.0 - normAlpha); // если normAlpha близка к 1, может быть большим

		do
		{
			// Вычисляем x_{k+1} = β + α * x_k
			for (int i = 0; i < n; i++)
			{
				double sum = 0.0;
				for (int j = 0; j < n; j++)
					sum += alpha[i, j] * xPrev[j];
				xCurr[i] = beta[i] + sum;
			}

			// Вычисляем норму разности
			double[] diff = new double[n];
			for (int i = 0; i < n; i++)
				diff[i] = xCurr[i] - xPrev[i];
			double diffNorm = Norm(diff, normType);

			// Проверка условия остановки (2.8)
			double tol = (normAlpha < 1.0) ? ((1.0 - normAlpha) / normAlpha * epsilon) : epsilon;
			if (diffNorm <= tol)
			{
				return xCurr;
			}

			// Подготовка к следующей итерации
			Array.Copy(xCurr, xPrev, n);
			iteration++;

			if (iteration > maxIterations)
				throw new Exception($"Превышено максимальное число итераций ({maxIterations}). Последняя норма разности: {diffNorm:F6}");

		} while (true);
	}

	private static (double[,] alpha, double[] beta) BuildAlphaBeta(double[,] A, double[] b)
	{
		int n = A.GetLength(0);
		double[,] alpha = new double[n, n];
		double[] beta = new double[n];

		for (int i = 0; i < n; i++)
		{
			double diag = A[i, i];
			if (Math.Abs(diag) < 1e-15)
				throw new InvalidOperationException($"Нулевой диагональный элемент в строке {i + 1}.");

			beta[i] = b[i] / diag;
			for (int j = 0; j < n; j++)
			{
				if (i == j)
					alpha[i, j] = 0.0;
				else
					alpha[i, j] = -A[i, j] / diag;
			}
		}
		return (alpha, beta);
	}

	/// <summary>
	/// Вычисляет норму вектора.
	/// </summary>
	public static double Norm(double[] vector, NormType type)
	{
		double result = 0;
		switch (type)
		{
			case NormType.L1:
				foreach (var v in vector) result += Math.Abs(v);
				break;
			case NormType.L2:
				foreach (var v in vector) result += v * v;
				result = Math.Sqrt(result);
				break;
			case NormType.LInf:
				result = double.NegativeInfinity;
				foreach (var v in vector) result = Math.Max(result, Math.Abs(v));
				break;
		}
		return result;
	}

	/// Вычисляет норму матрицы, согласованную с выбранной векторной нормой:
	/// L1 — максимум суммы модулей по столбцам,
	/// L2 — норма Фробениуса (корень из суммы квадратов всех элементов),
	/// LInf — максимум суммы модулей по строкам.
	public static double Norm(double[,] matrix, NormType type)
	{
		int n = matrix.GetLength(0);
		double result = 0;

		switch (type)
		{
			case NormType.L1:
				// максимум по столбцам
				for (int j = 0; j < n; j++)
				{
					double sum = 0;
					for (int i = 0; i < n; i++)
						sum += Math.Abs(matrix[i, j]);
					if (sum > result) result = sum;
				}
				break;

			case NormType.L2:
				// Фробениус
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
						result += matrix[i, j] * matrix[i, j];
				result = Math.Sqrt(result);
				break;

			case NormType.LInf:
				// максимум по строкам
				for (int i = 0; i < n; i++)
				{
					double sum = 0;
					for (int j = 0; j < n; j++)
						sum += Math.Abs(matrix[i, j]);
					if (sum > result) result = sum;
				}
				break;
		}
		return result;
	}

	/// Проверяет достаточное условие сходимости метода простых итераций (||α|| < 1).
	public static bool CheckConvergence(double[,] A, NormType normType = NormType.LInf)
	{
		int n = A.GetLength(0);
		double[,] alpha = new double[n, n];
		for (int i = 0; i < n; i++)
		{
			double diag = A[i, i];
			if (Math.Abs(diag) < 1e-15) return false; // нулевой диагональный элемент
			for (int j = 0; j < n; j++)
				alpha[i, j] = (i == j) ? 0.0 : -A[i, j] / diag;
		}
		return Norm(alpha, normType) < 1.0;
	}
}