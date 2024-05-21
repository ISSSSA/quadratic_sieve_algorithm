from random import randint
from .utils import gcd, product, isqrt, kth_iroot


def is_probable_prime(a):
    """ Оценка вероятности что число простое методом Миллера-Рабина
    """
    if a == 2:
        return True

    if a == 1 or a % 2 == 0:
        return False

    return rabin_miller_primality_test(a, 50)


def rabin_miller_primality_test(a, iterations):
    """ Тест простоты Рабин-Миллера
    """
    r, s = 0, a - 1

    while s % 2 == 0:
        r += 1
        s //= 2

    for _ in range(iterations):
        n = randint(2, a - 1)
        x = pow(n, s, a)
        if x == 1 or x == a - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, a)
            if x == a - 1:
                break
        else:
            return False
    return True


def check_perfect_power(n):
    """ Проверьте, является ли данное целое число совершенной степенью
    """
    prime = small_primes[-1]
    for p in small_primes:
        pth_root = kth_iroot(n, p)
        if pth_root < prime:
            break
        if pth_root ** p == n:
            return (pth_root, p)
    return None


def check_factor(n, i, factors):
    while n % i == 0:
        n //= i
        factors.append(i)
        if is_probable_prime(n):
            factors.append(n)
            n = 1
    return n


def find_small_primes(n, upper_bound):
    print("Попытка деления и инициализации маленьких простых чисел...")
    global small_primes
    is_prime = [True] * (upper_bound + 1)
    is_prime[0:2] = [False] * 2
    factors = []
    small_primes = []
    max_i = isqrt(upper_bound)
    rem = n
    for i in range(2, max_i + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

            for j in range(i ** 2, upper_bound + 1, i):
                is_prime[j] = False

    for i in range(max_i + 1, upper_bound + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

    print("Простые числа инициализированы!")
    return factors, rem


def find_prime_factors(n):
    print("Проверяем является ли {} идеальной силой ...".format(n))
    perfect_power = check_perfect_power(n)
    if perfect_power:
        print("{} это {}^{}".format(n, perfect_power[0], perfect_power[1]))
        factors = perfect_power[0]
    else:
        print("Not a perfect power")
        digits = len(str(n))
        if digits <= 30:
            print("Использую факторизацию методом р0 Полларда " + \
                  "для факторизации {} ({} чисел)".format(n, digits))
            factors = [brent_factorise(n)]
        else:
            print("Использование квадратичного решета для факторизации " + \
                  "{} ({} чисел)".format(n, digits))
            # factors = siqs_factorise(n)

    prime_factors = []
    for f in set(factors):
        for pf in find_all_prime_factors(f):
            prime_factors.append(pf)

    return prime_factors


def find_all_prime_factors(n):
    rem = n
    factors = []

    while rem > 1:
        if is_probable_prime(rem):
            factors.append(rem)
            break

        for f in find_prime_factors(rem):
            print("Найдены простые множители: {}".format(f))
            assert is_probable_prime(f)
            assert rem % f == 0
            while rem % f == 0:
                rem //= f
                factors.append(f)

    return factors


def _pollard_brent_func(c, n, x):
    y = (x ** 2) % n + c
    if y >= n:
        y -= n

    assert y >= 0 and y < n
    return y


def brent_factorise(n, iterations=None):
    y, c, m = (randint(1, n - 1) for _ in range(3))
    r, q, g = 1, 1, 1
    i = 0
    while g == 1:
        x = y
        for _ in range(r):
            y = _pollard_brent_func(c, n, y)
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = _pollard_brent_func(c, n, y)
                q = (q * abs(x - y)) % n
            g = gcd(q, n)
            k += m
        r *= 2
        if iterations:
            i += 1
            if i == iterations:
                return None

    if g == n:
        while True:
            ys = _pollard_brent_func(c, n, ys)
            g = gcd(abs(x - ys), n)
            if g > 1:
                break
    return g


def pollard_brent_iterator(n, factors):
    rem = n
    while True:
        if is_probable_prime(n):
            factors.append(n)
            rem = 1
            break

        digits = len(str(n))
        if digits < 45:
            iterations = 20
        else:
            iterations = 25

            f = brent_factorise(rem, iterations)
            if f and f < rem:
                if is_probable_prime(f):
                    print("Используя алгоритма метод-р0 Полларда для поиска простых множителей: " + \
                          "{}".format(f))
                    factors.append(f)
                    rem //= f
                else:
                    print("Brent's (Pollard's rho): Composite factor " + \
                          "found: {}".format(f))
                    rem_f = pollard_brent_iterator(f, factors)
                    rem = (rem // f) * rem_f
            else:
                print("Не найдено малых протсых множителей")
                break
    return rem


def factorise(n):
    if type(n) != int or n < 1:
        raise ValueError("Числа должны быть положительными")

    print("Факторизация {} ({} чисел)...".format(n, len(str(n))))

    if n == 1:
        return []

    if is_probable_prime(n):
        return [n]

    factors, rem = find_small_primes(n, 1000000)

    if factors:
        print("Prime factors found so far:")
        factors_temp = []
        for _ in factors:
            if _ not in factors_temp:
                factors_temp.append(_)
        print(*factors_temp, sep=', ')
    else:
        print("Не найдено маленьких простых множителей!")

    if rem != 1:
        digits = len(str(rem))
        if digits > 30:
            print("Используем р0-метод Полларда " + \
                  "для поиска большиз множителей...")
            rem = pollard_brent_iterator(rem, factors)
        if rem > 1:
            for f in find_all_prime_factors(rem):
                factors.append(f)

    factors.sort()
    assert product(factors) == n
    for p in factors:
        assert is_probable_prime(p)
    return factors


def main():
    n = int(input("Введите число для факторизации: "))
    result = factorise(n)
    new_result = []

    for _ in result:
        if _ not in new_result:
            new_result.append(_)

    print("\nПростые множители: {}".format(new_result))


if __name__ == '__main__':
    main()
