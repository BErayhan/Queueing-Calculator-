#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void MainMenu();
void MMCKK();
void MM1K();
void MM1();
void MMC();
double factorial(double);
void ending();

int main() { MainMenu(); }

void MainMenu() {
  system("cls");
  int choice;

  printf("Pilih Operasi queueing yang diinginkan:\n");
  printf("1.M/M/1 (single server)\n");
  printf("2.M/M/c (multi server, tapi severnya punya rentang jumlah)\n");
  printf("3.M/M/c/K/K (multi server dan arrival org nya eksponensial (k))\n");
  printf("4.M/M/1/k (single server, tapi arrivalnya dibatasi eksponen)\n");
  printf("\nPilihan Operasi : ");
  scanf("%d", &choice);

  switch (choice) {
  case 1:
    MM1();
    break;

  case 2:
    MMC();
    break;

  case 3:
    MMCKK();
    break;

  case 4:
    MM1K();
    break;

  default:
    printf("Pilihan Anda tidak tersedia");
  }

  return;
}

void MM1(void) {

  system("cls");
  system("color 1F");
  printf("~ KALKULATOR M/M/1 ~");

  // Variabel yang digunakan
  float lambda, mu, P0, L, Lq, W, Wq, U, I;

  // Meminta Input
  printf("\n\njumlah orang yang bisa dilayani setiap jam = ");
  scanf("%f", &mu);

  printf("jumlah orang datang setiap jam = ");
  scanf("%f", &lambda);

  // Rumus-rumus
  U = lambda / mu;
  P0 = 1 - (lambda / mu);
  L = lambda / (mu - lambda);
  Lq = L * U;
  W = 1 / (mu - lambda);
  Wq = W * U;
  I = 1 - U;

  // Error Handling
  if (P0 <= 0) {
    printf("\n!!! antrian ini telah memenuhi atau bahkan melebihi kapasitas "
           "server, menjadi rekomendasi untuk penambahan server !!!\n");
    P0 = 0;
    L = Lq = mu;
    W = Wq = 1;
    U = P0;
    I = 1 - U;
  }

  // Tampilan Hasil Perhitungan
  printf("\nProbabilitas tidak ada orang di dalam antrian (P0) = %.1f", P0);
  printf("\nRata-rata jumlah orang di seluruh sistem antrian setiap jam (L) = "
         "%.1f orang",
         L);
  printf(
      "\nRata-rata jumlah orang di garis antrian setiap jam (Lq) = %.1f orang",
      Lq);
  printf("\nRata-rata waktu tunggu setiap orang di seluruh sistem antrian (W) "
         "= %.1f jam",
         W);
  printf(
      "\nRata-rata waktu tunggu setiap orang di garis antrian (Wq) = %.1f jam",
      Wq);
  printf("\nProbabilitas Server sibuk/ada antrian (U) = %.1f", U);
  printf("\nProbabilitas Server kosong/tidak ada antrian (I) = %.1f", I);
  printf("\n\n");
  ending();
}

void MMCKK(void) {
  system("cls");
  system("color 3F");
  double lambda, mu, channel, K;
  double P[1000];
  int i;
  printf("~ KALKULATOR M/M/c/K/K ~\n\n");

  printf("Masukan Nilai Berikut\n");
  printf("LAMBDA (jumlah orang datang setiap jam) = ");
  scanf("%lf", &lambda);
  printf("\nMU (jumlah orang yang bisa dilayani setiap jam = ");
  scanf("%lf", &mu);
  printf("\nchannel: ");
  scanf("%lf", &channel);
  printf("\nk: ");
  scanf("%lf", &K);

  if (lambda < 0 || mu < 0) {
    printf("Lambda dan Mu harus Positif\nSilahkan kembali ke menu awal");
    getchar();
    MainMenu();
  } else if (K < channel) {
    printf("Jumlah K harus lebih besar sama dengan Channel \nSilahkan kembali "
           "ke menu awal");
    getchar();
    MainMenu();
  } else if (channel < 1) {
    printf("jumlah channel harus lebih besar sama dengan 1\nSilahkan kembali "
           "ke menu awal");
    getchar();
    MainMenu();
  }

  // Initialize variables
  double OfferedLoad = lambda / mu;
  double Kfac = factorial(K);
  double P0 = 1;

  // Calculate P0 and P(i)
  for (i = 1; i <= K; i++) {
    if (channel > 1 && i < channel) {
      P[i] = (Kfac / factorial(i) / factorial(K - i)) * pow(OfferedLoad, i);
    } else {
      P[i] = (Kfac / factorial(K - i) / factorial(channel) /
              pow(channel, i - channel)) *
             pow(OfferedLoad, i);
    }
    P0 += P[i];
  }

  P0 = 1 / P0;

  // Calculate other performance metrics
  double L = 0, LQ = 0, LambdaEffective = K * lambda * P0;
  for (i = 1; i <= K; i++) {
    P[i] *= P0;
    L += i * P[i];
    LQ += fmax(0, i - channel) * P[i];
    LambdaEffective += lambda * (K - i) * P[i];
  }

  double w = L / LambdaEffective;
  double wQ = LQ / LambdaEffective;
  double rho = LambdaEffective / channel / mu;

  // Output results
  printf("\n=== HASILNYA =====\n");
  printf("rho: %lf\n", rho);
  printf("L: %lf\n", L);
  printf("w: %lf\n", w);
  printf("wQ: %lf\n", wQ);
  printf("LQ: %lf\n", LQ);
  printf("P0: %lf\n", P0);
  printf("LambdaEffective: %lf\n\n", LambdaEffective);

  ending();
}

double factorial(double n) {
  // Implement factorial function
  double result = 1;
  int i;
  for (i = 2; i <= n; i++) {
    result *= i;
  }
  return result;
}

void ending() {
  int choice;
  printf(
      "\nApakah Anda ingin menggunakan kalkulator lagi? (1: Ya, 0: Tidak): ");
  scanf("%d", &choice);
  if (choice == 1) {
    MainMenu();
  }
}

void MMC() {
  // step kalkulasi disini

  float rho;
  float L;
  float W;
  float Wq;
  float Lq;
  float P0;

  system("cls");
  system("color 3F");
  double lambda, mu, channel, K;
  double P[1000];
  int i;

  printf("Masukan Nilai Berikut\n");
  printf("LAMBDA (jumlah orang datang setiap jam) = ");
  scanf("%lf", &lambda);
  printf("\nMU (jumlah orang yang bisa dilayani setiap jam = ");
  scanf("%lf", &mu);
  printf("\nchannel: ");
  scanf("%lf", &channel);

  if (lambda < 0 || mu < 0) {
    printf("Lambda dan Mu harus Positif\nSilahkan kembali ke menu awal");
    getchar();
    MainMenu();
  } else if (channel < 1) {
    printf("jumlah channel harus lebih besar sama dengan 1\nSilahkan kembali "
           "ke menu awal");
    getchar();
    MainMenu();
  }

  rho = lambda / mu / channel;
  float offered_load = lambda / mu;
  float factor = 1.;
  float p0 = 1.;

  for (i = 1; i < channel; i++) {
    factor *= offered_load / i;
    p0 += factor;
  }
  float cfactorial = tgamma(channel + 1.);
  p0 = p0 + factor * offered_load / channel / (1. - rho);
  p0 = 1. / p0;

  L = offered_load + pow(offered_load, channel + 1) * p0 / channel /
                         cfactorial / pow(1. - rho, 2);
  W = L / lambda;
  Wq = W - 1. / mu;
  Lq = Wq * lambda;

  printf("\n=== HASILNYA =====\n");
  printf("rho: %lf\n", rho);
  printf("L: %lf\n", L);
  printf("w: %lf\n", W);
  printf("wQ: %lf\n", Wq);
  printf("LQ: %lf\n", Lq);
  printf("P0: %lf\n", p0);
  ending();
}

void MM1K() {
  double Lambda, Mu, K, rho, OfferedLoad, Factor, P0, rhosum, cfactorial, PN, LQ, LambdaEffective, wQ, w, L;
	int i;
	double c = 1;
  // Input
  printf("Masukan Nilai Berikut\n");
  printf("LAMBDA (jumlah orang datang setiap jam) = ");
  scanf("%lf", &Lambda);
  printf("\nMU (jumlah orang yang bisa dilayani setiap jam = ");
  scanf("%lf", &Mu);
  printf("\nk: ");
  scanf("%lf", &K);

  // Validation
  if (Lambda < 0 || Mu < 0) {
    printf("Arrival rate and service rate must be positive\n");
    return;
  } else if (c < 1) {
    printf("Number of servers must be 1 or greater\n");
    return;
  } else if (K < c) {
    printf("Capacity must be at least as large as the number of servers\n");
    return;
  } else if (Lambda == c * Mu) {
    printf("This program does not handle the case Lambda equal c*Mu\n");
    return;
  }

  // Calculation
  rho = Lambda / Mu / c;
  OfferedLoad = Lambda / Mu;
  Factor = 1;
  P0 = 1;

  for (i = 1; i <= c; i++) {
    Factor *= OfferedLoad / i;
    P0 += Factor;
  }

  if (c < K) {
    rhosum = rho;
    if (c + 1 < K) {
      for (i = c + 2; i <= K; ++i) {
        rhosum += pow(rho, i - c);
      }
    }
    P0 += Factor * rhosum;
  }

  P0 = 1 / P0;
  cfactorial = tgamma(c + 1); // Using tgamma for factorial
  PN = pow(OfferedLoad, K) / cfactorial / pow(c, K - c) * P0;
  LQ = P0 * pow(OfferedLoad, c) * rho / cfactorial / pow(1 - rho, 2) *
       (1 - pow(rho, K - c) - (K - c) * pow(rho, K - c) * (1 - rho));
  LambdaEffective = Lambda * (1 - PN);
  wQ = LQ / LambdaEffective;
  w = wQ + 1 / Mu;
  L = LambdaEffective * w;

  // Output
  printf("\n=== HASILNYA =====\n");
  printf("rho: %lf\n", LambdaEffective / c / Mu);
  printf("L: %lf\n", L);
  printf("w: %lf\n", w);
  printf("wQ: %lf\n", wQ);
  printf("LQ: %lf\n", LQ);
  printf("P0: %lf\n", P0);
  printf("PN: %lf\n", PN);
  printf("LambdaEffective: %lf\n", LambdaEffective);
  ending();
}
