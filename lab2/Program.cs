namespace lab2;
// Вариант №5
class Programm 
{
	static void Main() 
	{
		// Матрица A
		var matrixA = new float[,]
		{
			{ 10.0f, 20.0f, 30.0f },
			{ 40.0f, 80.00001f, 60.0f },
			{ 5.0f, -15.0f, 25.0f }
		};

		// Вектор B
		var vectorB = new float[] { 60.0f, 180.00001f, 15.0f };

		Console.WriteLine("Исходная система:");
		Utility.PrintMatrix(matrixA);
		Utility.PrintVector(vectorB);

		float[] resNo = Gauss.GaussMethod(Utility.CopyMatrix(matrixA), Utility.CopyVector(vectorB));
		Console.WriteLine("Решение (x1, x2, x3):");
		Utility.PrintVector(resNo);

		float[] res = Gauss.GaussMethodWithKey(Utility.CopyMatrix(matrixA), Utility.CopyVector(vectorB));
		Console.WriteLine("Решение (x1, x2, x3):");
		Utility.PrintVector(res);

		// Матрица A
		var matrixA2 = new float[,]
		{
			{ 1f, 2f, 15f, 1f },
			{ -1f, 3f, 2f, 20f },
			{ 10f, 1f, 1f, 1f },
			{ 2f, 3f, 14f, -2f }
		};

		// Вектор B
		var vectorB2 = new float[] { 19f, 24f, 13f, 17f };

		SimpleIteration.Run(matrixA2, vectorB2);

		TestRandom(5);
	}

	static void TestRandom(int n)
	{
		Console.WriteLine($"\nСлучайный тест (n={n}) с диагональным преобладанием:");
		var rand = new Random();
		double[] a = new double[n - 1];
		double[] b = new double[n];
		double[] c = new double[n - 1];
		double[] f = new double[n];

		// Генерируем матрицу с диагональным преобладанием
		for (int i = 0; i < n; i++)
		{
			double sum = 0;
			if (i > 0) sum += Math.Abs(a[i - 1] = rand.NextDouble() * 0.5);
			if (i < n - 1) sum += Math.Abs(c[i] = rand.NextDouble() * 0.5);
			b[i] = sum + 1.0 + rand.NextDouble(); // диагональ больше суммы модулей внедиагональных
		}

		// Генерируем случайное решение и вычисляем правую часть
		double[] trueX = new double[n];
		for (int i = 0; i < n; i++) trueX[i] = rand.NextDouble() * 10 - 5;

		f[0] = b[0] * trueX[0] + c[0] * trueX[1];
		for (int i = 1; i < n - 1; i++)
			f[i] = a[i - 1] * trueX[i - 1] + b[i] * trueX[i] + c[i] * trueX[i + 1];
		f[n - 1] = a[n - 2] * trueX[n - 2] + b[n - 1] * trueX[n - 1];

		var computedX = TridiagonalSolver.Solve(a, b, c, f);

		Console.WriteLine("Заданное решение:  " + string.Join("  ", trueX));
		Console.WriteLine("Вычисленное решение: " + string.Join("  ", computedX));

		double maxError = 0;
		for (int i = 0; i < n; i++)
			maxError = Math.Max(maxError, Math.Abs(trueX[i] - computedX[i]));
		Console.WriteLine($"Максимальная ошибка: {maxError:e4}");
	}
}