#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 256

//離散フーリエ変換
void dft1d(float rx[], float ix[], float Ru[], float Iu[]) {
    for (int u = 0; u < N; u++) {
        Ru[u] = 0;
        Iu[u] = 0;
        for (int x = 0; x < N; x++) {
            Ru[u] += rx[x] * cos((2 * M_PI * u * x) / (float)N) + ix[x] * sin((2 * M_PI * u * x) / (float)N);
            Iu[u] += ix[x] * cos((2 * M_PI * u * x) / (float)N) - rx[x] * sin((2 * M_PI * u * x) / (float)N);
        }
    }
}

//左上と右下、右上と左下を入れ替える
void swap(float a[][N], float b[][N],float x[][N],float y[][N]) {
    for (int i = 0; i < N/2; i++) {
        for (int j = 0; j < N/2; j++) {
                x[i + N / 2][j + N / 2] = a[i][j];
                x[i][j + N / 2] = a[i + N / 2][j];
                x[i + N / 2][j] = a[i][j + N / 2];
                x[i][j] = a[i + N / 2][j + N / 2];

                y[i + N / 2][j + N / 2] = b[i][j];
                y[i][j + N / 2] = b[i + N / 2][j];
                y[i + N / 2][j] = b[i][j + N / 2];
                y[i][j] = b[i + N / 2][j + N / 2];
        }
    }
}

//パワースペクトル作成
void make_power_spectrum(float Pu[][N],float R[][N],float I[][N]) {
    for (int i = 0; i < N; i++) {       
        for (int j = 0; j < N; j++) {
            Pu[i][j] = R[i][j] * R[i][j] + I[i][j] * I[i][j];
            Pu[i][j] += 1;
        }
    }

    for (int i = 0; i < N; i++) {       //パワースペクトルを対数変換
        for (int j = 0; j < N; j++) {
            Pu[i][j] = log(Pu[i][j]);
        }
    }
}

//入力された半径の円の1/4を作成
void make_quarter(int f[][N / 2], int a) {
    for (int i = 0; i < N / 2; i++) {
        for (int j = 0; j < N / 2; j++) {
            float n = sqrt(i * i + j * j);
            if (n <= a) {
                f[i][j] = 1;
            }
            else {
                f[i][j] = 0;
            }
        }
    }
}

//1/4円から円を作成
void make_circle(int x[][N / 2], int f[][N]) {
    for (int i = 0; i < N / 2; i++) {
        for (int j = 0; j < N / 2; j++) {
            f[i][j] = x[N / 2 - i - 1][N / 2 - j - 1];
            f[i][j + N / 2] = x[N / 2 - i - 1][j];
            f[i + N / 2][j] = x[i][N / 2 - j - 1];
            f[i + N / 2][j + N / 2] = x[i][j];
        }
    }
}

int main()
{
    int D;
    unsigned char fr[N][N], fi[N][N];    //入力画像
    float gr[N][N], gi[N][N];           //ワーク画像
    float rxy[N], ixy[N];               //抽出配列
    float Ru[N], Iu[N];                 //出力配列
    float R[N][N], I[N][N];             //フーリエ結果,逆フーリエ結果
    float R1[N][N], I1[N][N];           //入れ替え後
    float R2[N][N], I2[N][N];           //フィルタ処理後　
    float Pu[N][N];                     //パワースペクトル出力用
    unsigned char fR[N][N], fI[N][N];   //結果出力用
    int quarter[N / 2][N / 2];          //1/4円生成用
    int F[N][N];                        //理想低域通過フィルタ
    FILE* fp1, * fp2, * fp3, * fp4;
    
    fp1 = fopen("noiseimage.raw", "rb");//画像読み込み
    fread(fr, 1, N * N, fp1);
    fclose(fp1);

    for (int i = 0; i < N; i++) {       //虚部画像作成
        for (int j = 0; j < N; j++) {
            fi[i][j] = 0;
        }
    }

    for (int i = 0; i < N; i++) {       //行方向の離散フーリエ変換
        for (int j = 0; j < N; j++) {
            rxy[j] = (float)fr[i][j];
            ixy[j] = (float)fi[i][j];
        }
        dft1d(rxy, ixy, Ru, Iu);
        for (int j = 0; j < N; j++) {
            gr[i][j] = Ru[j];
            gi[i][j] = Iu[j];
        }
    }

    for (int i = 0; i < N; i++) {       //列方向の離散フーリエ変換
        for (int j = 0; j < N; j++) {
            rxy[j] = gr[j][i];
            ixy[j] = gi[j][i];
        }
        dft1d(rxy, ixy, Ru, Iu);
        for (int j = 0; j < N; j++) {
            R[j][i] = Ru[j];
            I[j][i] = Iu[j];
        }
    }

    swap(R, I, R1, I1);                 //直流分を中央にする

    make_power_spectrum(Pu, R1, I1);    //パワースペクトル作成

    fp2 = fopen("power_spectrum.raw", "wb");//パワースペクトル出力
    fwrite(Pu, 4, N * N, fp2);
    fclose(fp2);
    
    printf("D=");
    scanf("%d", &D);

    make_quarter(quarter, D);
    make_circle(quarter, F);
    
    for (int i = 0; i < N; i++) {       //フィルタリング処理
        for (int j = 0; j < N; j++) {
            R2[i][j] = R1[i][j] * F[i][j];
            I2[i][j] = I1[i][j] * F[i][j];
        }
    }

    swap(R2, I2, R1, I1);               //中央にした直流分元の位置に戻す

    for (int i = 0; i < N; i++) {       //行方向の逆離散フーリエ変換
        for (int j = 0; j < N; j++) {
            rxy[j] = R1[i][j];
            ixy[j] = I1[i][j];
        }
        dft1d(ixy, rxy, Iu, Ru);
        for (int j = 0; j < N; j++) {
            gr[i][j] = Ru[j]/N;
            gi[i][j] = Iu[j]/N;
        }
    }

    for (int i = 0; i < N; i++) {       //列方向の逆離散フーリエ変換
        for (int j = 0; j < N; j++) {
            rxy[j] = gr[j][i];
            ixy[j] = gi[j][i];
        }
        dft1d(ixy, rxy, Iu, Ru);
        for (int j = 0; j < N; j++) {
            R[j][i] = Ru[j]/N;
            I[j][i] = Iu[j]/N;
        }
    }
    
    for (int i = 0; i < N; i++) {       //逆フーリエ変換の結果の絶対値を取り、負の値を除去
        for (int j = 0; j < N; j++) {
            R[i][j] = fabs(R[i][j]);
            I[i][j] = fabs(I[i][j]);
        }
    }
  
    for (int i = 0; i < N; i++) {       //出力のためunsigned charに変換
        for (int j = 0; j < N; j++) {
            fR[i][j] = (unsigned char)R[i][j];
            fI[i][j] = (unsigned char)I[i][j];
        }
    }
    //実画像を保存
    fp3 = fopen("R_image.raw", "wb");
    fwrite(fR, 1, N * N, fp3);
    fclose(fp3);
    //虚画像を保存
    fp4 = fopen("I_image.raw", "wb");
    fwrite(fI, 1, N * N, fp4);
    fclose(fp4);
    return 0;
}