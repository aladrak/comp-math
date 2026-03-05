namespace lab2;
static class Utility
{
    public static void PrintMatrix(float[,] m)
    {
        int n = m.GetLength(0);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                Console.Write($"{m[i, j],8:F3} ");
            Console.WriteLine();
        }
    }

	//static void PrintMatrix(float[,] a, float[] b)
	//{
	//	int n = b.Length;
	//	for (int i = 0; i < n; i++)
	//	{
	//		for (int j = 0; j < n; j++)
	//		{
	//			Console.Write($"{a[i, j],8:F5} * x{j + 1}");
	//			if (j < n - 1) Console.Write(" + ");
	//			else Console.Write($" = {b[i],8:F5}");
	//		}
	//		Console.WriteLine();
	//	}
	//}

	public static void PrintVector(float[] x)
    {
		Console.WriteLine("[ " + string.Join(", ", x) + " ]");
		//for (int i = 0; i < x.Length; i++)
		//	Console.WriteLine($"x{i + 1} = {x[i]}");
	}

	public static float[,] CopyMatrix(float[,] original)
	{
		int rows = original.GetLength(0);
		int cols = original.GetLength(1);
		float[,] copy = new float[rows, cols];
		Array.Copy(original, copy, original.Length);
		return copy;
	}

	public static float[] CopyVector(float[] original) => (float[])original.Clone();
}