
def product(lst):

    prod = 1
    for _ in lst:
        prod *= _
    return prod


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def kth_iroot(n, k):
    u = n
    s = n + 1
    while u < s:
        s = u
        t = (k - 1) * s + n // pow(s, k - 1)
        u = t // k
    return s


def isqrt(n):
    if n < 0:
        raise ValueError("Квадратынй корень отрицательного числа!")
    x = int(n)
    if n == 0:
        return 0
    a, b = divmod(x.bit_length(), 2)
    n = 2 ** (a + b)
    while True:
        y = (n + x // n) // 2
        if y >= x:
            return x
        x = y