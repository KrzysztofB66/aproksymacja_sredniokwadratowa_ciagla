import math


# Dane wejściowe
def funkcja(x):
    return math.sqrt(3 + 2*math.pow(x, 2))


a = -1
b = 1
m = 2  # stopień wielomianu
x = -0.25


def eliminacja_gausa(A, b):
    n = len(A)
    M = A

    i = 0
    for x in M:
        x.append(b[i])
        i += 1

    for k in range(n):
        for i in range(k, n):
            if abs(M[i][k]) > abs(M[k][k]):
                M[k], M[i] = M[i], M[k]
            else:
                pass

        for j in range(k+1, n):
            q = float(M[j][k]) / M[k][k]
            for m in range(k, n+1):
                M[j][m] -= q * M[k][m]

    x = [0 for i in range(n)]

    x[n-1] = round(float(M[n-1][n])/M[n-1][n-1], 6)
    for i in range(n-1, -1, -1):
        z = 0
        for j in range(i+1, n):
            z = z + float(M[i][j])*x[j]
        x[i] = round(float(M[i][n] - z)/M[i][i], 6)
    return x


def metoda_simpsona(f, a, b, n):
    h = (b - a) / n
    x0 = f(a) + f(b)
    x1 = sum(f(a + i * h) for i in range(1, n, 2))
    x2 = sum(f(a + i * h) for i in range(2, n-1, 2))

    return (h / 3) * (x0 + 4*x1 + 2*x2)


def funkcje_bazowe(m):
    return [lambda x, i=i: x**i for i in range(m+1)]


def aproksymacja_sredniokwadratowa(f, a, b, m, xx):
    phi = funkcje_bazowe(m)

   # Obliczanie macierzy A i wektora b
    m_A = []
    w_b = []
    for i in range(m+1):
        m_A.append([])
        w_b.append(metoda_simpsona(lambda y: funkcja(y) * phi[i](y), a, b, 20))
        for j in range(m+1):
            m_A[i].append(metoda_simpsona(lambda x: phi[i](x) * phi[j](x), a, b, 20))

    a_n = eliminacja_gausa(m_A, w_b)
    wynik = 0
    for i in range(m+1):
        wynik += a_n[i]*phi[i](x)

    return wynik


print(aproksymacja_sredniokwadratowa(funkcja, a, b, m, x))







