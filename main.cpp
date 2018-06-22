#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;
void heapSort(int *&priority, int *&perm, int &n);

// algorytm szeregowania zadan metodą wstawień; wersja bazowa z akceleracją
// zwraca cmax uszeregowania (nie globalnie optymalne)
long algorytmNEH(int &n, int &m, int **&p, long &cmax, int *&perm);

int main() {
    int n, m;
    long cmax;
    int **p;
    int *perm;
    int reps=1;
    chrono::duration<double> czas_obliczen_NEH = chrono::duration<double>::zero();
    chrono::duration<double> elapsed_seconds;

    string plik_wyj = "permutacje.txt";
    string plik_wej = "neh.data.txt";
    ifstream plik;
    ofstream wyjscie;
    string linia;
    wyjscie.open(plik_wyj, ofstream::trunc);
    wyjscie.close();
    plik.open(plik_wej);

    while (!plik.eof()) {
        getline(plik, linia);
        if (linia[0] == 'd' && linia[1] == 'a') { // napotkano linie 'data', czytaj od nastepnej
            plik >> n >> m; // wczytaj rozmiar
            // zarezerwuj pamiec
            p = new int *[m + 1];
            for (int i = 1; i <= m; ++i)
                p[i] = new int[n + 1];
            // wczytaj pozostale dane
            for (int j = 1; j <= n; ++j)
                for (int i = 1; i <= m; ++i)
                    plik >> p[i][j];

            perm = new int[n + 1];
            chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            for (int i=0; i<reps; ++i) // licz czas wykonania dla 'reps' powtorzen algorytmu
                algorytmNEH(n, m, p, cmax, perm);
            chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            elapsed_seconds = end - start;
            czas_obliczen_NEH += elapsed_seconds;

            wyjscie.open(plik_wyj, ofstream::app);
            //for (int j = 1; j <= n; ++j) wyjscie << perm[j] << " ";
            wyjscie << endl << cmax ;//<< endl;
            //wyjscie << n << " " << m << " "; // rozmiar danych, jesli potrzebne
            wyjscie << /*"Czas obliczen: " <<*/ elapsed_seconds.count() << endl;
            wyjscie.close();

            for (int i = 1; i <= m; ++i)
                delete[] p[i];
            delete[] p;
            delete[] perm;
        }
    }
    plik.close();
    wyjscie.open(plik_wyj, ofstream::app);
    //wyjscie << "Calkowity czas obliczen: " << czas_obliczen_NEH.count() << "s" << endl;
    cout << "Calkowity czas obliczen: " << czas_obliczen_NEH.count() << "s" << endl;
    wyjscie.close();
}

long algorytmNEH(int &n, int &m, int **&p, long &cmax, int *&perm) {
    int *priority = new int[n + 1]; // priorytet sortowania
    int i, j, l, t;                 // zmienne pomocnicze
    int **r = new int *[m + 1];     // odleglosc do wierzcholka (1,1) grafu
    for (i = 0; i <= m; ++i) r[i] = new int[n + 1];
    int **q = new int *[m + 2];     // dlugosci drogi do wierzch (m,i) grafu
    for (i = 0; i <= m + 1; ++i) q[i] = new int[n + 2];
    int *dist = new int[m + 1];
    long c; // biezaca wartosc funkcji celu

    for (i = 1; i <= n; ++i) {    //liczenie priorytetow
        t = 0;
        for (j = 1; j <= m; ++j)
            t += p[j][i];
        priority[i] = t * n - i;  // suma czasow wykonania wszystkich zadan pomniejszona o indeks
    }
    for (i = 1; i <= n; ++i)     // stworz permutacje bazowa, zeby moc sortowac
        perm[i] = i;

    heapSort(priority, perm, n); //sortowanie przez kopcowanie wzgledem priorytetu

    for (i = 0; i <= m; ++i) r[i][0] = 0;  // zadanie 0 nie istnieje
    for (j = 0; j <= n; ++j) r[0][j] = 0;  // maszyna 0 nie istnieje
    for (j = 0; j <= n + 1; ++j) q[m + 1][j] = 0; // maszyna m+1 nie istnieje
    dist[0] = 0; // ze wzgledu na brak definicji

    t = 1;  // pozycja kolejnego zadania (na poczatku 1)
    for (int zad = 2; zad <= n; ++zad) {
        // odswiez r i q
        for (i = 1; i <= m; ++i)
            for (j = t; j <= zad; ++j)  // do zad, zamiast zad-1; aby szybicej policzyć cmax dla kolejnego zadania na koncu
                r[i][j] = max(r[i - 1][j], r[i][j - 1]) + p[i][perm[j]]; // licz odleglosci do (1,1)

        for (i = 0; i <= m; ++i) q[i][zad] = 0; //prawy brzeg 'q' przed dodaniem zadania 'zad'

        for (i = m; i >= 1; --i)
            for (j = zad - 1; j >= t; --j)  // przepisz q z j-1 na j od pozycji t
                q[i][j] = q[i][j - 1];

        for (i = m; i >= 1; --i)            // odswiez q od pozycji t do 1
            for (j = t; j >= 1; --j)
                q[i][j] = max(q[i][j + 1], q[i + 1][j]) + p[i][perm[j]];

        // wyznacz najlepsza pozycje w permutacji
        cmax = r[m][zad];                           // cmax najdluzsza sciezka w grafie
        t = zad;                                    // pozycja kolejnego zadania
        for (j = zad - 1; j >= 1; --j) {    // wyznacz najlepsza pozycje dla nowego zadania i wylicz cmax
            for (l = 1; l <= m; ++l)        // wyznacz D(l), droge z po wstawieniu zadania w j'te miejsce
                dist[l] = max(dist[l - 1], r[l][j - 1]) + p[l][perm[zad]];
            c = dist[1] + q[1][j];
            for (l = 2; l <= m; ++l)
                c = max(c, (long) dist[l] + q[l][j]); // oblicz maksymalna dlugosc uszeregowania dla permutacji
            if (c <= cmax) {    // jesli maksymalna dlugosc krotsza od dotychczasowego cmax, zapamietaj
                cmax = c;
                t = j;  // pozycja dla kolejnego zadania , jesli nie znajdziesz w petli, tzn, ze pozycja na koncu permutacji
            }
        }
        // zapisz permutacje
        i = perm[zad];
        for (j = zad; j > t; --j)
            perm[j] = perm[j - 1];
        perm[t] = i;
    }
    // zwolnienie pamieci
    for (i = 0; i <= m + 1; ++i)
        delete[] q[i];
    delete[] q;
    for (i = 0; i <= m; ++i)
        delete[] r[i];
    delete[] r;
    delete[] dist;
    delete[] priority;
    return cmax;
}


void heapSort(int *&priority, int *&perm, int &n) {
    int p_ind;  // rodzic (parent index)
    int ch_ind; //  dziecko (child index)
    int gch_ind; // wieksze dziecko (greater child index)
    int x, y;

    for (int i = 2; i <= n; i++) { // budowanie kopca
        ch_ind = i;
        p_ind = ch_ind / 2;
        x = priority[i];
        y = perm[i];

        while ((p_ind > 0) && (priority[p_ind] > x)) {
            perm[ch_ind] = perm[p_ind];
            priority[ch_ind] = priority[p_ind];
            ch_ind = p_ind;
            p_ind = (ch_ind >> 1);
        }
        priority[ch_ind] = x;
        perm[ch_ind] = y;
    }

    for (int i = n; i > 1; i--)      // rozebranie kopca
    {
        swap(priority[1], priority[i]); // element najwiekszy na koniec
        swap(perm[1], perm[i]);
        p_ind = 1;
        ch_ind = 2;
        while (ch_ind < i) //z
        {
            if ((ch_ind + 1 < i) &&
                (priority[ch_ind + 1] < priority[ch_ind])) //jesli istnieje prawe dziecko i jest mniejsze od lewego
                gch_ind = ch_ind + 1; // mniejsze prawe
            else
                gch_ind = ch_ind;   // mniejsze lewe
            if (priority[gch_ind] >= priority[p_ind])
                break; // jesli mniejsze dziecko jest wieksze od rodzica przerwij 'while'
            swap(priority[p_ind], priority[gch_ind]); // jesli nie, zamien rodzica z wiekszym dzieckiem
            swap(perm[p_ind], perm[gch_ind]);
            p_ind = gch_ind;
            ch_ind = (p_ind << 1);
        }
    } // zmniejsz rozmiar kopca
}
