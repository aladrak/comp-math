namespace lab2;
static class Gauss
{
	public static float[] GaussMethod(float[,] a, float[] b)
	{
		Console.WriteLine("\nМетод Гаусса БЕЗ выбора главного элемента");
		int n = b.Length;
		for (int k = 0; k < n - 1; k++)
		{
			if (Math.Abs(a[k, k]) < 1e-9f)
			{
				Console.WriteLine("Ошибка: Деление на ноль");
				continue;
			}

			for (int i = k + 1; i < n; i++)
			{
				float factor = a[i, k] / a[k, k];
				for (int j = k; j < n; j++) a[i, j] -= factor * a[k, j];
				b[i] -= factor * b[k];
			}
		}
		// Обратный ход
		float[] x = new float[n];
		for (int i = n - 1; i >= 0; i--)
		{
			float sum = 0;
			for (int j = i + 1; j < n; j++) sum += a[i, j] * x[j];
			x[i] = (b[i] - sum) / a[i, i];
		}
		return x;
	}

	public static float[] GaussMethodWithKey(float[,] a, float[] b)
	{
		Console.WriteLine("\nМетод Гаусса с выбором по строке");
		int n = b.Length;
		for (int k = 0; k < n - 1; k++)
		{
			int maxRow = k;
			float maxVal = Math.Abs(a[k, k]);
			// Ищем строку с максимальным по модулю элементом в текущем столбце k
			for (int i = k + 1; i < n; i++)
			{
				if (Math.Abs(a[i, k]) > maxVal)
				{
					maxVal = Math.Abs(a[i, k]);
					maxRow = i;
				}
			}
			// Меняем местами текущую строку k и строку maxRow
			if (maxRow != k)
			{
				// Меняем местами строки матрицы
				for (int j = 0; j < n; j++)
				{
					float temp = a[k, j];
					a[k, j] = a[maxRow, j];
					a[maxRow, j] = temp;
				}
				// Меняем местами элементы вектора b
				float tempB = b[k];
				b[k] = b[maxRow];
				b[maxRow] = tempB;
			}
			// Исключение неизвестных
			for (int i = k + 1; i < n; i++)
			{
				float factor = a[i, k] / a[k, k];
				for (int j = k; j < n; j++) a[i, j] -= factor * a[k, j];
				b[i] -= factor * b[k];
			}
		}
		// Обратный ход
		float[] x = new float[n];
		for (int i = n - 1; i >= 0; i--)
		{
			float sum = 0;
			for (int j = i + 1; j < n; j++) sum += a[i, j] * x[j];
			x[i] = (b[i] - sum) / a[i, i];
		}
		return x;
	}
}