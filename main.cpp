#include <iostream>
#include <vector>
#include <math.h>


std::vector<double> gauss(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();
    
    // Tworzenie rozszerzonej macierzy
    for (int i = 0; i < n; i++)
        A[i].push_back(b[i]);
    
    // Eliminacja
    for (int i = 0; i < n; i++) {
        // Szukanie maksimum w kolumnie
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
            if (abs(A[k][i]) > abs(A[maxRow][i]))
                maxRow = k;
        
        swap(A[i], A[maxRow]);

        // Normalizacja wiersza
        double diag = A[i][i];
        for (int j = i; j <= n; j++)
            A[i][j] /= diag;
        
        // Eliminacja dolna
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i];
            for (int j = i; j <= n; j++)
                A[k][j] -= factor * A[i][j];
        }
    }
    
    // Podstawienie wsteczne
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
    }
    return x;
}

// Funkcja aproksymująca wielomianem 6. stopnia
std::vector<double> least_squares_polynomial_6(const std::vector<double>& xs, const std::vector<double>& ys) {
    int n = xs.size();
    // if (n != 7) {
    //     throw std::runtime_error("Potrzeba dokładnie 7 punktów!");
    // }

    const int degree = 6;
    std::vector<std::vector<double>> A(degree + 1, std::vector<double>(degree + 1, 0.0));
    std::vector<double> b(degree + 1, 0.0);

    // Wypełnianie macierzy normalnej
    for (int i = 0; i <= degree; i++) {
        for (int j = 0; j <= degree; j++) {
            for (int k = 0; k < n; k++)
                A[i][j] += pow(xs[k], i + j);
        }
        for (int k = 0; k < n; k++)
            b[i] += ys[k] * pow(xs[k], i);
    }

    // Rozwiązywanie układu
    return gauss(A, b);
}

double calc_polynomial(std::vector<double> coeffs, double x){
    double sum = 0;
    for(int i = 0; i < coeffs.size(); i++){
        sum += pow(x, i) * coeffs[i];
    }

    return sum;
}

int main(int argc, char const *argv[])
{
   /*  // h0 - poczatkowa dlugosc materialu
    const double h = 1;
    
    // ### 3
    double start_stress = 0;

    // ### 4 i 5
    int best_index = find_closest_point_stress(start_stress);
    double stress_ps = y_stress[best_index];
    double strain_ps = x_strain[best_index];
    double absolute_deformation_ps = h * strain_ps;    // dl_ps = h0 * epsilon[%]

    // 6
    double l0 = h + absolute_deformation_ps; // l0 = h + dl_ps

    // 7 wykres epsilon(naprezenia)

    // 8    E = 10%, 30%, 60%
    double EpsilonL_1 = 0.1;
    double EpsilonL_2 = 0.3;
    double EpsilonL_3 = 0.6;

    // dl = E_l * l0 + dl_ps
    // wartości x dla danych odkształceń El (epsilon l)
    double new_dl_1 = EpsilonL_1 * l0 + absolute_deformation_ps;
    double new_dl_2 = EpsilonL_2 * l0 + absolute_deformation_ps;
    double new_dl_3 = EpsilonL_3 * l0 + absolute_deformation_ps;

    // delta sigma - wartość y z wykresu dla kolejnych El (epsilon l)
    double sigma_1 = y_stress[find_closest_point_strain(new_dl_1)];
    double sigma_2 = y_stress[find_closest_point_strain(new_dl_2)];
    double sigma_3 = y_stress[find_closest_point_strain(new_dl_3)];


    // MODULY SIECZNE
    double dSigma1 = sigma_1 - stress_ps;
    double dSigma2 = sigma_2 - stress_ps;
    double dSigma3 = sigma_3 - stress_ps;
    double E_1 = dSigma1 / EpsilonL_1;
    double E_1 = dSigma2 / EpsilonL_2;
    double E_1 = dSigma3 / EpsilonL_3; */

    // odksztalcenie  strain    [%]
    std::vector<double> x_values = { 
        10.558997642412209, 20.263961250036015, 30.745340190475133,
        40.450303798098936, 50.93164127443481, 60.24841794770258, 70.34160995378562
    };

    // naprezenie - stress      [MPa]
    std::vector<double> y_values = {
        9.833505087516746, 9.139774153357092, 8.29958890172689,
        7.806269117176516, 7.536484956482864, 7.467111945398491, 7.582733699148773
    };


    const double h = 1;

    // wyliczenie wartości początkowych dla zadanego odkształcenia wstępnego = 10%
    double naprezenie_pocz      = y_values[0] * 0.01;  // zamiana z % na liczbe
    double odksztalcenie_pocz   = x_values[0] * 0.01;  // zamiana z % na liczbe

    double l_0 = h + h*odksztalcenie_pocz;      // h + dl

    // korekta dla zadanego odksztalcenie poczatkowe
    std::vector<double> x_corrected;
    std::vector<double> y_corrected;
    x_corrected.reserve( x_values.size());
    y_corrected.reserve( x_values.size());

    for(int i = 0; i < x_values.size(); i++){
        y_corrected.push_back(y_values[i] - naprezenie_pocz);
        x_corrected.push_back(x_values[i] - odksztalcenie_pocz);
    }

    // wyznaczanie modułów siecznych
    // moduly sieczne   20%, 30%, 40%, 50%, 60%, 70%
    std::vector<double> EpsilonL = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    std::vector<double> El;

    for(int i = 0; i < EpsilonL.size(); i++){
        // dSigma = sigma - sigma_ps - co zostało policzone wczesniej -> y_corrected
        El.push_back(y_corrected[i+1]/EpsilonL[i]); // dSigma / EpsilonL
    }

    // a0 = coeffs[0] = wyraz wolny
    std::vector<double> coeffs = least_squares_polynomial_6(x_corrected, y_corrected);

    printf("\nmaprox coeffs:\n");
    for(int i = 0; i < coeffs.size(); i++){
        printf("a%d: %0.15lf, ", i,coeffs[i]);
    }
    // a6 + a5 + a4 + a3 + a2 + a1 + a0 
    // ^6 + ^5 + ^4 + ^3 + ^2 + ^1 + ^0 
    // pochodna
    // 6*a6 + 5*a5
    // ^5 + ^4 + ^3 + ^2 + ^1 + ^0 

    // wyliczenie wspolczynnikow wielopmianu pochodnej bazujac na fukncji aproksymujacej
    std::vector<double> coeffs_der = std::vector(6, 0.0);

    printf("\nderivative coeffs:\n");
    for(int i = 0; i < coeffs_der.size(); i++){
        coeffs_der[i] = i * coeffs[i];
        printf("a%d: %0.15lf", i, coeffs_der[i]);
    }

    // wyznaczenie modulow stycznych

    std::vector<double> E_st;
    // podstawienie do pochodnej w punkcie EpsilonL
    for(int i = 1; i < x_corrected.size(); i++){
        E_st.push_back(calc_polynomial(coeffs_der, x_corrected[i]));
    }



    // wypisanie wyników
    printf("\n### MODULY SIECZNE:\n");
    for(int i = 0; i < El.size(); i++){
        // EpsilonL
        printf("%2.0lf\%: E = %lf\n", EpsilonL[i]*100, El[i]);

    }
    printf("\n### MODULY STYCZNE:\n");
    for(int i = 0; i < E_st.size(); i++){
        // EpsilonL
        printf("%2.0lf\%: E_st = %lf\n", EpsilonL[i]*100, E_st[i]);

    }

    return 0;
}
