/*
 * Лабораторная работа №2
 * Solution систем линейных алгебраических уравнений
 * Вариант 5
 * Задание 1. Метод Гаусса без выбора главного элемента и с выбором главного элемента ПО СТРОКАМ (вариант 5) — для SLAE 1.
 * Задание 2. Метод Простых Iter — для SLAE 2.
 * Задание 3. Метод прогонки для произвольной трёхдиагональной системы.
 * g++ -O2 -std=c++17 -o lab2_var5 lab2_var5.cpp
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

typedef vector<float>         VecF;
typedef vector<vector<float>> MatF;

void printMatrix(const MatF& A, const VecF& b, const string& title = "")
{
    if (!title.empty())
        cout << "\n" << title << "\n";
    int n = (int)A.size();
    cout << fixed << setprecision(6);
    for (int i = 0; i < n; i++) {
        cout << "  |";
        for (int j = 0; j < n; j++)
            cout << setw(12) << A[i][j];
        cout << " |" << setw(12) << b[i] << " |\n";
    }
}

void printVector(const VecF& x, const string& title = "")
{
    if (!title.empty())
        cout << "\n" << title << "\n";
    cout << fixed << setprecision(6);
    for (int i = 0; i < (int)x.size(); i++)
        cout << "  x[" << i + 1 << "] = " << x[i] << "\n";
}

float normInf(const VecF& v)
{
    float m = 0.f;
    for (float vi : v) m = max(m, fabsf(vi));
    return m;
}

float normDiff(const VecF& a, const VecF& b)
{
    VecF d(a.size());
    for (int i = 0; i < (int)a.size(); i++) d[i] = a[i] - b[i];
    return normInf(d);
}

float matNormInf(const MatF& M)
{
    int n = (int)M.size();
    float mx = 0.f;
    for (int i = 0; i < n; i++) {
        float s = 0.f;
        for (int j = 0; j < n; j++) s += fabsf(M[i][j]);
        mx = max(mx, s);
    }
    return mx;
}

void extractTridiagonal(const MatF& A, VecF& c, VecF& d, VecF& e) {
    int n = (int)A.size();
    c.assign(n, 0.f);
    d.assign(n, 0.f);
    e.assign(n, 0.f);
    for (int i = 0; i < n; i++) {
        d[i] = A[i][i];
        if (i > 0) c[i] = A[i][i - 1];
        if (i < n - 1) e[i] = A[i][i + 1];
    }
}

// Метод Гаусса БЕЗ выбора главного элемента
VecF gaussNoSelect(MatF A, VecF b)
{
    int n = (int)A.size();
    for (int k = 0; k < n - 1; k++) {
        if (fabsf(A[k][k]) < 1e-9f)
            cout << "  [WARNING] Small leading element in step k=" << k
                 << " (" << A[k][k] << ")\n";
        for (int m = k + 1; m < n; m++) {
            float factor = A[m][k] / A[k][k];
            for (int l = k; l < n; l++) A[m][l] -= factor * A[k][l];
            b[m] -= factor * b[k];
        }
    }
    VecF x(n);
    for (int i = n - 1; i >= 0; i--) {
        float s = b[i];
        for (int j = i + 1; j < n; j++) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return x;
}

// Метод Гаусса с выбором главного элемента ПО СТРОКЕ (вариант 5)
//   На k-м шаге ищем max|A[k][j]|, j=k..n-1; переставляем столбцы.
VecF gaussRowSelect(MatF A, VecF b)
{
    int n = (int)A.size();
    vector<int> col(n);
    for (int i = 0; i < n; i++) col[i] = i;

    for (int k = 0; k < n - 1; k++) {
        int mc = k;
        float mv = fabsf(A[k][k]);
        for (int j = k + 1; j < n; j++)
            if (fabsf(A[k][j]) > mv) { mv = fabsf(A[k][j]); mc = j; }
        if (mc != k) {
            for (int i = 0; i < n; i++) swap(A[i][k], A[i][mc]);
            swap(col[k], col[mc]);
            cout << "  k=" << k << ": swap col " << k << "<->" << mc
                 << "  pivot=" << mv << "\n";
        }
        for (int m = k + 1; m < n; m++) {
            float f = A[m][k] / A[k][k];
            for (int l = k; l < n; l++) A[m][l] -= f * A[k][l];
            b[m] -= f * b[k];
        }
    }
    VecF xp(n);
    for (int i = n - 1; i >= 0; i--) {
        float s = b[i];
        for (int j = i + 1; j < n; j++) s -= A[i][j] * xp[j];
        xp[i] = s / A[i][i];
    }
    VecF x(n);
    for (int i = 0; i < n; i++) x[col[i]] = xp[i];
    return x;
}

// Построение итерационной формы x = Beta*x + gamma
void buildIterForm(const MatF& A, const VecF& b, MatF& Beta, VecF& gamma)
{
    int n = (int)A.size();
    Beta.assign(n, VecF(n, 0.f));
    gamma.assign(n, 0.f);
    for (int i = 0; i < n; i++) {
        gamma[i] = b[i] / A[i][i];
        for (int j = 0; j < n; j++)
            if (j != i) Beta[i][j] = -A[i][j] / A[i][i];
    }
}

// Функция для проверки сходимости метода простых итераций
void checkConvergence(const MatF& A, const VecF& b) {
    int n = (int)A.size();
    cout << "\n--- Convergence analysis ---\n";
    
    // Проверка диагонального преобладания
    bool diagDominant = true;
    for (int i = 0; i < n; i++) {
        float diag = fabsf(A[i][i]);
        float sum = 0.f;
        for (int j = 0; j < n; j++) if (j != i) sum += fabsf(A[i][j]);
        if (diag <= sum) diagDominant = false;
        cout << "  Row " << i+1 << ": |" << diag << "| " 
             << (diag > sum ? ">" : "<=") << " " << sum 
             << (diag > sum ? " -> OK" : " -> NOT dominant") << "\n";
    }
    if (diagDominant) 
        cout << "  Matrix has diagonal dominance (sufficient condition for some iterative methods).\n";
    else 
        cout << "  Matrix does NOT have diagonal dominance.\n";
    
    // Построение Beta и вычисление её нормы
    MatF Beta;
    VecF gamma;
    buildIterForm(A, b, Beta, gamma);
    float normBeta = matNormInf(Beta);
    cout << "  ||Beta||_inf = " << fixed << setprecision(6) << normBeta << "\n";
    if (normBeta < 1.f)
        cout << "  Sufficient condition for convergence is satisfied.\n";
    else
        cout << "  [WARNING] ||Beta|| >= 1, method may diverge.\n";
}

// Функция для перестановки строк с целью улучшения диагонального преобладания
void improveDiagonalDominance(MatF& A, VecF& b) {
    int n = (int)A.size();
    vector<int> perm(n);
    for (int i = 0; i < n; i++) perm[i] = i;
    
    // Жадный алгоритм: для каждой позиции i ищем строку с максимальным диагональным элементом среди оставшихся
    for (int i = 0; i < n; i++) {
        int best = i;
        float maxDiag = fabsf(A[perm[i]][i]);
        for (int j = i+1; j < n; j++) {
            if (fabsf(A[perm[j]][i]) > maxDiag) {
                maxDiag = fabsf(A[perm[j]][i]);
                best = j;
            }
        }
        if (best != i) swap(perm[i], perm[best]);
    }
    
    // Применяем перестановку строк
    MatF newA(n, VecF(n));
    VecF newb(n);
    for (int i = 0; i < n; i++) {
        newA[i] = A[perm[i]];
        newb[i] = b[perm[i]];
    }
    A = move(newA);
    b = move(newb);
    
    cout << "  Rows reordered to improve diagonal dominance.\n";
}

// Обновлённая функция simpleIteration с возможностью автоматического преобразования
VecF simpleIteration(const MatF& A, const VecF& b, float eps, int& iters,
                     bool verbose = true, bool autoTransform = true)
{
    MatF curA = A;
    VecF curb = b;
    int n = (int)curA.size();
    
    // Автоматическое преобразование, если необходимо
    if (autoTransform) {
        MatF Beta; VecF gamma;
        buildIterForm(curA, curb, Beta, gamma);
        float normBeta = matNormInf(Beta);
        if (normBeta >= 1.f) {
            cout << "  Attempting to improve convergence by reordering rows...\n";
            improveDiagonalDominance(curA, curb);
        }
    }
    
    // Построение итерационной формы x = Beta*x + gamma
    MatF Beta(n, VecF(n, 0.f));
    VecF gamma(n, 0.f);
    
    for (int i = 0; i < n; i++) {
        if (fabsf(curA[i][i]) < 1e-12f) {
            cout << "  [ERROR] Diagonal element A[" << i << "][" << i 
                 << "] = " << curA[i][i] << " close to zero!\n";
            iters = -1;
            return VecF();
        }
        gamma[i] = curb[i] / curA[i][i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                Beta[i][j] = -curA[i][j] / curA[i][i];
            }
        }
    }
    
    float normBeta = matNormInf(Beta);
    if (verbose) {
        cout << fixed << setprecision(6);
        cout << "  Matrix norm Beta (inf-norm): " << normBeta << "\n";
        if (normBeta < 1.f) {
            cout << "  The sufficient condition for convergence is satisfied (||Beta|| < 1).\n";
        } else {
            cout << "  [WARNING] Norm Beta >= 1, the method may not converge.\n";
        }
        cout << "\n  k  |  ||x^{k+1} - x^{k}||_inf\n"
                "  ---+------------------------\n";
    }
    
    // Начальное приближение x^{(0)} = gamma
    VecF x = gamma;
    VecF xOld(n);
    iters = 0;
    
    for (int iter = 0; iter < 10000; iter++) {
        xOld = x;
        for (int i = 0; i < n; i++) {
            float sum = 0.f;
            for (int j = 0; j < n; j++) {
                sum += Beta[i][j] * xOld[j];
            }
            x[i] = gamma[i] + sum;
        }
        iters++;
        float diff = normDiff(x, xOld);
        if (verbose) {
            cout << "  " << setw(3) << iters << " |  " << diff << "\n";
        }
        if (diff < eps) {
            if (verbose) cout << "\n  Convergence achieved in " << iters << " Iter.\n";
            break;
        }
        if (!isfinite(diff) || diff > 1e10f) {
            if (verbose) cout << "  [ERROR] The method diverges!\n";
            break;
        }
    }
    return x;
}

// Метод прогонки
//   c[i]*x[i-1] + d[i]*x[i] + e[i]*x[i+1] = f[i]
//   c[0]=0, e[n-1]=0
bool isTridiagonal(const MatF& A, float eps = 1e-6f) {
    int n = (int)A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (abs(i - j) > 1 && fabsf(A[i][j]) > eps) {
                return false;
            }
        }
    }
    return true;
}

// Основная функция метода прогонки, принимающая диагонали
VecF thomasAlgorithm(const VecF& c, const VecF& d, const VecF& e, const VecF& f) {
    int n = (int)d.size();
    VecF alpha(n), beta(n);
    
    if (fabsf(d[0]) < 1e-12f) {
        cout << "  [ERROR] d[0] = " << d[0] << " is zero!\n";
        return VecF();
    }
    alpha[0] = -e[0] / d[0];
    beta[0]  =  f[0] / d[0];
    
    for (int i = 1; i < n; ++i) {
        float denom = d[i] + c[i] * alpha[i-1];
        if (fabsf(denom) < 1e-12f) {
            cout << "  [ERROR] denom[" << i << "] = " << denom << " is zero!\n";
            return VecF();
        }
        alpha[i] = (i < n-1) ? -e[i] / denom : 0.f;
        beta[i]  = (f[i] - c[i] * beta[i-1]) / denom;
    }
    
    VecF x(n);
    x[n-1] = beta[n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = alpha[i] * x[i+1] + beta[i];
    }
    return x;
}

VecF thomasSolver(const MatF& A, const VecF& b) {
    int n = (int)A.size();
    
    // Проверка, что матрица квадратная и размеры совпадают
    if (n == 0 || A[0].size() != n || b.size() != n) {
        cout << "  [ERROR] Invalid matrix/vector dimensions.\n";
        return VecF();
    }
    
    if (!isTridiagonal(A)) {
        cout << "  [WARNING] The matrix is NOT tridiagonal!\n";
    }
    
    // Извлечение диагоналей
    VecF c(n, 0.f), d(n, 0.f), e(n, 0.f);
    for (int i = 0; i < n; ++i) {
        d[i] = A[i][i];
        if (i > 0) c[i] = A[i][i-1];
        if (i < n-1) e[i] = A[i][i+1];
    }
    
    return thomasAlgorithm(c, d, e, b);
}

int main()
{
    cout << " LabWork No.2  |  Var 5\n"
         << " Solution SLAE: Gauss, SimpleIteration, Run-through\n";

    // SLAE 1
    cout << "\n=== SLAE 1 ===\n";
    MatF A1 = {
        {20.f,  2.f,       15.f, 1.f },
        {-1.f,  25.f,  2.f, 20.f},
        { 10.f,  1.f,       15.f, 1.f },
        {2., 3., 14., 20.}
    };
    VecF b1 = {19.f, 24.f, 13.f, 17.};

    printMatrix(A1, b1, "Initial system:");

    cout << "\n--- Gaussian method WITHOUT selecting the Main element ---\n";
    VecF x1n = gaussNoSelect(A1, b1);
    printVector(x1n, "Solution:");

    cout << "\n--- Gaussian method with selection BY ROW ---\n";
    VecF x1r = gaussRowSelect(A1, b1);
    printVector(x1r, "Solution:");

    cout << "\n  Difference in solutions ||x_row - x_nosel||_inf = "
         << fixed << setprecision(6) << normDiff(x1r, x1n) << "\n";

/*
    Пояснение: строки 1 и 2 исходной матрицы почти пропорциональны
    (строка2 ~ 4*строка1), поэтому на шаге k=1 ведущий элемент
    после первого исключения ~1e-5. Деление на малое число
    усиливает ошибки округления 32-битной арифметики.
    Выбор главного элемента по строке переставляет столбцы,
    помещая наибольший элемент на позицию ведущего, что
    повышает численную устойчивость.
*/
    
    // SLAE 2 (МЕТОД ПРОСТЫХ ИТЕРАЦИЙ)
    cout << "\n=== SLAE 2 (Simple Iteraion) ===\n";
    MatF A2orig = {
        {20.f,  2.f,       15.f, 1.f },
        {-1.f,  25.f,  2.f, 20.f},
        { 10.f,  1.f,       15.f, 1.f },
        {2., 3., 14., 20.}
    };
    VecF b2orig = {19.f, 24.f, 13.f, 17.};

    printMatrix(A2orig, b2orig, "Initial system:");

    MatF A2 = {A2orig[2], A2orig[3], A2orig[0], A2orig[1]};
    VecF b2 = {b2orig[2], b2orig[3], b2orig[0], b2orig[1]};

    printMatrix(A2, b2, "Transformed system (lines in order 3,4,1,2):");

    cout << "\n  Diagonal dominance of the transformed:\n";
    for (int i = 0; i < 4; i++) {
        float diag = fabsf(A2[i][i]), rest = 0.f;
        for (int j = 0; j < 4; j++) if (j != i) rest += fabsf(A2[i][j]);
        cout << "    i=" << i+1 << ": |" << diag << "| vs "
             << rest << (diag > rest ? "  OK" : "  FAIL") << "\n";
    }

    cout << "\n--- Simple iteration method (eps=1e-3) ---\n";
    int iters2 = 0;
    VecF x2 = simpleIteration(A2, b2, 1e-3f, iters2);
    printVector(x2, "Result:");
    cout << "  Iter: " << iters2 << "\n";

    cout << "\n--- Check (Gauss with selection) ---\n";
    VecF x2ref = gaussRowSelect(A2orig, b2orig);
    printVector(x2ref, "Reference:");
    cout << "  Error: " << fixed << setprecision(8)
         << normDiff(x2, x2ref) << "\n";

    cout << "\n--- Selecting the minimum error ---\n";
    for (float ep : {1e-1f, 1e-2f, 1e-3f, 1e-4f, 1e-5f}) {
        int it = 0;
        VecF xt = simpleIteration(A2, b2, ep, it, false);
        cout << fixed << setprecision(4)
             << "  eps=" << ep << "  iters=" << setw(3) << it
             << "  err=" << setprecision(8) << normDiff(xt, x2ref) << "\n";
    }

    // Метод прогонки
    cout << "\n=== Run-through method (ex 3) ===\n";
    {
        MatF A = {
            {20.f, 2.f, 15.f, 1.f },
            {-1.f, 25.f,  2.f, 20.f},
            {10.f, 1.f, 15.f, 1.f },
            {2.f,  3.f, 14.f, 20.f}
        };
        VecF b = {19.f, 24.f, 13.f, 17.f};
        
        cout << "Original system:\n";
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                cout << setw(8) << A[i][j] << " ";
            cout << "| " << b[i] << "\n";
        }
        
        // Решение методом прогонки (через обёртку)
        VecF x = thomasSolver(A, b);
        if (!x.empty()) {
            printVector(x, "Solution using Thomas algorithm (extracted tridiagonal):");
        }
    }
}