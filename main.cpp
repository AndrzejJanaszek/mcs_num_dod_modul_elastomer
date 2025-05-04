#include <iostream>
#include <vector>
#include <math.h>

std::vector<double> gauss(std::vector<std::vector<double>> A, std::vector<double> b)
{
    int n = A.size();

    // Tworzenie rozszerzonej macierzy A|B
    for (int i = 0; i < n; i++)
        A[i].push_back(b[i]);

    // dla każdej kolumny
    for (int i = 0; i < n; i++)
    {
        // max w kolumnie i podmiana
        int max_row = i;
        for (int k = i + 1; k < n; k++)
            if (abs(A[k][i]) > abs(A[max_row][i]))
                max_row = k;

        swap(A[i], A[max_row]);

        // Normalizacja wiersza - pierwszy element 1 i reszta el. podzielona
        // przez jego początkową wartość
        // 1 a1/a0 a2/a0 ...
        // for el in row    // od diagonali
        double diag = A[i][i];
        for (int j = i; j <= n; j++)
            A[i][j] /= diag;

        // Eliminacja -> wszystko pod elementem "aktywnym" na zero
        // reszta przeliczona:
        // i-ty wiersz mnożymy razy element pod "aktywnym" i odejmujemy o k-tego wiersza
        for (int k = i + 1; k < n; k++)
        {
            double factor = A[k][i];
            for (int j = i; j <= n; j++)
                A[k][j] -= factor * A[i][j];
        }
    }

    // Podstawienie wsteczne   
    //                      itd             itd
    // 1*xn_2 + ax_n-1 + a*xn | b     => x_n-2 = b - a*xn - a*x_n-1
    //          1x_n-1 + a*xn | b     => x_n-1 = b - a*xn
    //                   1*xn | b     => xn = b
    // od "dołu" układu wyliczamy kolejne niewadome
    // zaczynając od tej którą możemy wyliczyć bezpośrednio
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
    }
    return x;
}

// Funkcja aproksymująca wielomianem 6. stopnia
std::vector<double> least_squares_polynomial_6(const std::vector<double> &xs, const std::vector<double> &ys)
{
    int n = xs.size();

    const int degree = 6;
    std::vector<std::vector<double>> A(degree + 1, std::vector<double>(degree + 1, 0.0));
    std::vector<double> b(degree + 1, 0.0);

    // Wypełnianie macierzy normalnej
    // macierz - układ równań powstaje z:
    // S_i = [ f(x) - W_6 ] ^2
    // gdzie W_6 -> a_0 + a_1*x + a_2*x^2 .... a_6*x^6
    // mamy 7 punktów więc będzie 7 równań liniowych czyli docelowo macierz 7x8 bo A|B
    for (int i = 0; i <= degree; i++)
    {
        for (int j = 0; j <= degree; j++)
        {
            for (int k = 0; k < n; k++)
                A[i][j] += pow(xs[k], i + j);
        }
        for (int k = 0; k < n; k++)
            b[i] += ys[k] * pow(xs[k], i);
    }

    // Rozwiązywanie układu
    return gauss(A, b);
}

double calc_polynomial(std::vector<double> coeffs, double x)
{
    double sum = 0;
    for (int i = 0; i < coeffs.size(); i++)
    {
        sum += pow(x, i) * coeffs[i];
    }

    return sum;
}

int main(int argc, char const *argv[])
{
    // *uwaga* punkt (0,0) - jest brany dla odkształcenia wstępnego 10%

    // odksztalcenie - strain    [%]
    std::vector<double> x_values = {
        10.643999745028083, 20.304102923051225, 30.679778357867743, 40.33988153589087, 40.33988153589087, 51.073327813018125, 60.37567543481654, 70.75135086963306};

    // naprezenie - stress      [MPa]
    std::vector<double> y_values = {
        13.270139553050981, 27.203794815976156, 43.98104414388616, 53.93365070698346, 53.93365070698346, 59.19431460994837, 60.616114679735276, 58.34123274569954};

    const double h = 1;

    // wyliczenie wartości początkowych dla zadanego odkształcenia wstępnego = 10%
    double naprezenie_pocz = y_values[0];
    double odksztalcenie_pocz = x_values[0];

    //                                    zmiana jednostki % na liczbę
    double l_0 = h + h * odksztalcenie_pocz * 0.01; // h + dl

    // korekta dla zadanego odksztalcenia poczatkowego
    // przesunięcie wykresu do punktu (0,0)
    std::vector<double> x_corrected;
    std::vector<double> y_corrected;
    x_corrected.reserve(x_values.size());
    y_corrected.reserve(x_values.size());

    for (int i = 0; i < x_values.size(); i++)
    {
        y_corrected.push_back(y_values[i] - naprezenie_pocz);
        x_corrected.push_back(x_values[i] - odksztalcenie_pocz);
    }

    printf("### Wartości po przeliczeniu:\n");
    printf("x:\n");
    for (int i = 0; i < x_corrected.size(); i++)
    {
        printf("%0.15lf, ", x_corrected[i]);
    }
    printf("\ny:\n");
    for (int i = 0; i < y_corrected.size(); i++)
    {
        printf("%0.15lf, ", y_corrected[i]);
    }

    // wyznaczanie modułów siecznych
    // moduly sieczne   20%, 30%, 40%, 50%, 60%, 70%
    std::vector<double> EpsilonL = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    std::vector<double> El;

    for (int i = 0; i < EpsilonL.size(); i++)
    {
        // dSigma = sigma - sigma_ps - co zostało policzone wczesniej -> y_corrected
        El.push_back(y_corrected[i + 1] / EpsilonL[i]); // dSigma / EpsilonL
    }

    // a0 = coeffs[0] = wyraz wolny
    std::vector<double> coeffs = least_squares_polynomial_6(x_corrected, y_corrected);

    printf("\n\n### Współczynniki wielomianu aproksymującego oraz pochodnej:\n");
    printf("maprox coeffs:\n");
    printf("a0-a6:\n");
    for (int i = 0; i < coeffs.size(); i++)
    {
        printf("%0.15lf, ", coeffs[i]);
    }
    // a6 + a5 + a4 + a3 + a2 + a1 + a0
    // ^6 + ^5 + ^4 + ^3 + ^2 + ^1 + ^0
    // pochodna
    // 6*a6 + 5*a5
    // ^5 + ^4 + ^3 + ^2 + ^1 + ^0

    // wyliczenie wspolczynnikow wielopmianu pochodnej bazujac na fukncji aproksymujacej
    std::vector<double> coeffs_der = std::vector(6, 0.0);

    printf("\nderivative coeffs:\n");
    printf("a0-a6:\n");
    for (int i = 0; i < coeffs_der.size(); i++)
    {
        coeffs_der[i] = i * coeffs[i];
        printf("%0.15lf", coeffs_der[i]);
    }

    // wyznaczenie modulow stycznych

    std::vector<double> E_st;
    // podstawienie do pochodnej w punkcie EpsilonL
    for (int i = 1; i < x_corrected.size(); i++)
    {
        E_st.push_back(calc_polynomial(coeffs_der, x_corrected[i]));
    }

    // wypisanie wyników
    printf("\n\n### MODULY SIECZNE:\n");
    for (int i = 0; i < El.size(); i++)
    {
        // EpsilonL
        printf("%2.0lf\%: E = %lf\n", EpsilonL[i] * 100, El[i]);
    }
    printf("\n### MODULY STYCZNE:\n");
    for (int i = 0; i < E_st.size(); i++)
    {
        // EpsilonL
        printf("%2.0lf\%: E_st = %lf\n", EpsilonL[i] * 100, E_st[i]);
    }

    return 0;
}
