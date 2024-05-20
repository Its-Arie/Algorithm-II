#include <iostream>
#include <complex>
using namespace std;
double pi = acos(-1);
typedef complex<double> cd;
cd *fft(int n, cd *s, bool invert)
{
    cd *a = new cd[n];
    if (n <= 1)
    {
        a[0] = s[0];
        return a;
    }
    cd *x = new cd[n / 2];
    cd *y = new cd[n / 2];
    for (int i = 0; i < n / 2; i++)
    {
        x[i] = s[2 * i];
        y[i] = s[2 * i + 1];
    }
    cd *b = fft(n / 2, x, invert);
    cd *c = fft(n / 2, y, invert);
    double power = -2 * pi / n * (invert ? -1 : 1);
    cd w(cos(power), sin(power));
    cd p = 1.0;
    for (int k = 0; k < n / 2; k++)
    {
        cd temp = p * c[k];
        a[k] = b[k] + temp;
        a[k + n / 2] = b[k] - temp;
        if (invert)
        {
            a[k] /= 2;
            a[k + n / 2] /= 2;
        }
        p = p * w;
    }
    return a;
}
int main()
{
    int n1, n2;
    cout << "size of polynomial 1 : ";
    cin >> n1;
    cd *p1 = new cd[n1];
    for (int i = 0; i < n1; i++)
    {
        double data;
        cout << "input for x^" << i << " : ";
        cin >> data;
        p1[i] = {data, 0};
    }
    cout << "size of polynomial 2 : ";
    cin >> n2;
    cd *p2 = new cd[n2];
    for (int i = 0; i < n2; i++)
    {
        double data;
        cout << "input for x^" << i << " : ";
        cin >> data;
        p2[i] = {data, 0};
    }
    int n = 1;
    while (n < n1 + n2 - 1)
    {
        n = n << 1;
    }
    cd *a = new cd[n];
    cd *b = new cd[n];
    for (int i = 0; i < n; i++)
    {
        if (i < n1)
        {
            a[i] = p1[i];
        }
        else
        {
            a[i] = 0;
        }
    }
    for (int i = 0; i < n; i++)
    {
        if (i < n2)
        {
            b[i] = p2[i];
        }
        else
        {
            b[i] = 0;
        }
    }
    cd *fa = fft(n, a, false);
    cd *fb = fft(n, b, false);
    for (int i = 0; i < n; i++)
    {
        fa[i] *= fb[i];
    }
    cd *ifft = fft(n, fa, true);
    cout << "the result is\n";
    for (int i = 0; i < n1 + n2 - 1; ++i)
    {
        cout << "coefficient of x^" << i << ": " << ifft[i].real() << endl;
    }
}